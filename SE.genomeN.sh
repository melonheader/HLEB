#!/usr/bin/env bash
### INFO: Count bases N per genomic feature for single end reads
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
			echo "-n, --n_cores           number of cores to parallelize jobs (default 6)"
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
				n_cores=6
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
results_path=$cnvs_path"/genomeN"

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
echo "Purpose: Count genome N in .bam files from Single-End Slamseq run"
echo "${#path_to_bams[@]} .bam files will be processed"
echo "Output dir: $results_path"



#{MAIN}
# define main function $1 - bam_path, $2 - base to count
countb () {
	s_time="$(date -u +%s)"
	# get bam basename and sample_name
	bam_file="$(basename -- $1)"
	sample_name="${bam_file%_S*}"
	echo ""
	echo "Processing sample $sample_name....."
	echo "Counting $2 background bases....."
	if [[ ! -d "./$2" ]]; then
		mkdir "./$2"
	fi
	cd $2
	
	# +
	samtools view -f 16 -b -h -L $bed_path $1 | \
	samtools mpileup -A -B -Q 1 -d 1000000 -f $fasta_path - | \
	awk -v base="$2" -v lbase="${bloop[$2]}" '($3==lbase || $3==base){print $1,$2,$3,$4,$5,"-"}' OFS="\t" > \
	cov"$2"."$sample_name".pileup
	# -
	samtools view -F 16 -b -h -L $bed_path $1 | \
	samtools mpileup -A -B -Q 1 -d 1000000 -f $fasta_path - | \
	awk -v base="$2" -v lbase="${bloop[$2]}" '($3==lbase || $3==base){print $1,$2,$3,$4,$5,"-"}' OFS="\t" >> \
	cov"$2"."$sample_name".pileup
	
	# filter unique
	awk '{gsub(">","<",$5); print}' cov"$2"."$sample_name".pileup | \
	awk '{print $1"_"$2"_"$3"_"$6,$4-gsub(/</,"",$5)}' OFS="\t" | \
	awk '($2>0)' | \
	awk '{a[$1]+=$2;}END{for(i in a)print i"\t"a[i];}' | \
	sed 's/_/\t/g' | \
	awk '{print $1,$2,$2+1,$3,$5,$4}' OFS="\t" > cov"$2"."$sample_name".bed

	# sort
	sort -k1,1 -k2,2n cov"$2"."$sample_name".bed > sorted.cov"$2"."$sample_name".bed

	# quantify conversions per gene
	intersectBed -s -wa -wb -sorted -a $bed_path -b sorted.cov"$2"."$sample_name".bed > cov"$2".gfeat."$sample_name".bed
	awk '{print $7,$8,$12,$13,$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$7}' OFS="\t" cov"$2".gfeat."$sample_name".bed | \
	awk '!seen[$5]++' |  \
	awk '{a[$1]+=$4;}END{for(i in a)print i"\t"a[i];}' > counts"$2".perGene."$sample_name".txt

	# report
	e_time="$(date -u +%s)"
	el_s=$(($e_time - $s_time))
	el_time=$(date --date='@'$el_s +%H:%M:%S)
	echo "genome background $2 in sample $sample_name processed"
	echo "Time elapsed: $el_time"
}


# declase a hash of bpairs
declare -A bloop=( ["T"]="t" ["C"]="c" ["G"]="g" ["A"]="a" )
# start
# set total start time
ts_time="$(date -u +%s)"
# relocate to the output directory
cd $results_path
# loop over samples
for x in ${path_to_bams[@]}; do
	# loop over background bases
	for y in ${!bloop[@]}; do
		## for timepoints later than 0 only do T
		timepoint_test=$(grep -q T0 $x)
		if [ -z $timepoint_test ] && [ $y != "T" ]; then
			continue
		fi
		(
			countb $x $y
		) &
		# allow parallel execution of N jobs
		if [[ $(jobs -r -p | wc -l) -gt $((n_cores - 1))  ]]; then
			# wait for a batch to finish
			wait
		fi
	done
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