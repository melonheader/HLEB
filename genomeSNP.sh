#!/usr/bin/env bash
### INFO: Call SNPs from provided bam files
### DATE: 02.08.2022
### AUTHOR: Artem Baranovskii, Miha Milek


#{SETUP}
set -e
# set variables
# ------------------------------ #
## Global
### Process arguments
while (( "$#" )); do
    case "$1" in
        -h|--help)
            echo "options:"
            echo "-h, --help              show brief help"
            echo "-i, --path_to_bams      path(s) to input bam files"
            echo "-o, --out_path          specify a directory to write results"
			echo "-b, --bed_path 		  path to bed with genomic features"
			echo "-f, --fasta_path        path to genome fasta corresponding to the bed file provided"
			echo "-n, --n_cores           number of cores to parallelize jobs (default 4)"
            exit 0
        ;;
        -i|--path_to_bams)
            shift 
            if test $# -gt 0; then
                path_to_bams=()
				args=( "$@" )
				set -- "${args[@]}"
				while (( $# )); do
					if [ ${1:0:1} == "-" ]; then
						break
					fi
					#echo "Path: $1"
					path_to_bams+=($1)
					shift
				done
				unset args
            else
                echo "No input bam files provided"
				exit 1
            fi
        ;;
        -o|--out_path)
            shift
            if test $# -gt 0; then
                out_path=$1
            else
                echo "No output dir specified"
                exit 1
            fi
            shift
		;;
		-b|--bed_path)
			shift 
			if test $# -gt 0; then
				bed_path=$1
			else
				echo "No path to bed file with genomic features provided"
				exit 1
			fi
			shift
        ;;
		-f|--fasta_path)
			shift 
			if test $# -gt 0; then
				fasta_path=$1
			else
				echo "No path to genomic fasta provided"
				exit 1
			fi
			shift
        ;;
		-v|--varscan_path)
			shift
			if test $# -gt 0; then
				varscan_path=$1
			else
				echo "No path to varscan executable provided"
				exit 1
			fi
			shift
		;;
		# -e|--experiment_name)
		# 	shift 
		# 	if test $# -gt 0; then
		# 		experiment_name=$1
		# 	else 
		# 		experiment_name="$(date -u +'%d.%m.%Y')_Run"
		# 	fi
		# 	shift
		# ;;
		-n|--n_cores)
			shift
			if test $# -gt 0; then
				n_cores=$1
			else
				n_cores=4
			fi
			shift
		;;
        *)
            echo "bad option"
            exit 1
        ;;
    esac
done



#{CONTROL}
# ----------------------------------------------------------------- #
## Paths
### source location
source=${BASH_source[0]}
while [ -L "$source" ]; do # resolve $source until the file is no longer a symlink
  source_dir=$( cd -P "$( dirname "$source" )" >/dev/null 2>&1 && pwd )
  source=$(readlink "$source")
  [[ $source != /* ]] && source=$source_dir/$source # if $source was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
source_dir=$( cd -P "$( dirname "$source" )" >/dev/null 2>&1 && pwd )

### Auxiliary
cntr_path=$source_dir"/counters"
cnvs_path=$out_path"/conversions"
results_path=$cnvs_path"/SNPs"

## Check for dirs to write results into
if [[ ! -d "$cnvs_path" ]]; then
	mkdir -p "$cnvs_path"
fi
if [[ ! -d "$results_path" ]]; then
	mkdir -p "$results_path"
fi
# check for .pl counter
if [[ ! -d "$cntr_path" ]]; then
	echo ""
	echo "$cntr_path"
	echo "^ .pl counters folder is not found. exiting....."
	exit 1
fi




# ----------------------------------------------------------------- #
# Echo settings
echo "Settings:"
echo "Purpose: Call SNPs from Slamseq run"
echo "${#path_to_bams[@]} .bam files will be processed"
echo ""



#{MAIN}
# define main function $1 - bam_path, $2 - base to count
call_snp () {
	s_time="$(date -u +%s)"
	# get bam basename and sample_name
	bam_file="$(basename -- $1)"
	sample_name="${bam_file%_trimmed*}"
	echo "Processing sample $sample_name....."

	samtools mpileup -B -A -q 255 -f $fasta_path $1 | \
	java -jar $varscan_path  mpileup2snp --strand-filter 0 --output-vcf --min-var-freq $minVarFreq --min-coverage $minCov --variants 1 | \
	awk '($4=="T"&&$5=="C")' | \
	awk '{print $1,$2,$2+1,$4$5,$10}' OFS="\t" | \
	sed 's/:/\t/g' | \
	awk '{print $1,$2,$3,$4,$10,"+"}' OFS="\t" > \
	counts."$sample_name".bed

	samtools mpileup -B -A -q 255 -f $fasta_path $1 | \
	java -jar $varscan_path  mpileup2snp --strand-filter 0 --output-vcf --min-var-freq $minVarFreq --min-coverage $minCov --variants 1 | \
	awk '($4=="A"&&$5=="G")' | \
	awk '{print $1,$2,$2+1,$4$5,$10}' OFS="\t" | \
	sed 's/:/\t/g' | \
	awk '{print $1,$2,$3,$4,$10,"-"}' OFS="\t" >> \
	counts."$sample_name".bed

	sort -k1,1 -k2,2n counts."$sample_name".bed > sorted.counts."$sample_name".bed
	intersectBed -s -wa -wb -sorted -a $bed_path -b sorted.counts."$sample_name".bed | \
	awk '($6=="+"&&$12=="TC")||($6=="-"&&$12=="AG")' > counts.gfeat."$sample_name".bed
	
	awk '{print $4,$8,$12,$13,$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$8}' OFS="\t" counts.gfeat."$sample_name".bed | \
	awk '!seen[$5]++' | \
	awk '{a[$1"_"$2]+=$4;}END{for(i in a)print i"\t"a[i];}' > countsTC.perGene."$sample_name".txt

	rm -f counts."$sample_name".bed

	# report
	e_time="$(date -u +%s)"
	el_s=$(($e_time - $s_time))
	el_time=$(date --date='@'$el_s +%H:%M:%S)
	echo "Finished SNP calling in $sample_name"
	echo "Time elapsed: $el_time"
}

# define threshold for varscan
minVarFreq="0.8"
minCov="10"
# start
# set total start time
ts_time="$(date -u +%s)"
# relocate to the output directory
cd $results_path
# loop over samples
for x in ${path_to_bams[@]}; do
	(
		call_snp $x
	) &
	# allow parallel execution of N jobs
	if [[ $(jobs -r -p | wc -l) -gt $((n_cores - 1))  ]]; then
		# wait for a batch to finish
		wait
	fi
done
# wait for unfinished jobs from the last batch
wait
# t report
te_time="$(date -u +%s)"
tel_s=$(($te_time - $ts_time))
tel_time=$(date --date='@'$tel_s +%H:%M:%S)
echo ""
echo "${#path_to_bams[@]} files procced"
echo "Total elapsed time: $tel_time"
echo "Done!"




