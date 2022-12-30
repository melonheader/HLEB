#!/usr/bin/env bash
### INFO: Count genome N -> N conversions from Slamseq bam files
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
results_path=$cnvs_path"/genomeNN"

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
echo "Experiment: $experiment_name"
echo "Purpose: Count genome N -> N conversions in .bam files from SPaired-End Slamseq run"
echo "${#path_to_bams[@]} .bam files will be processed"
echo ""



#{MAIN}
# define main function $1 - bam_path, $2 - main base, $3 - rev base
countc() {
	# set process start time
	s_time="$(date -u +%s)"
	# get bam basename and sample_name
	bam_file="$(basename -- $1)"
	sample_name="${bam_file%_trimmed*}"
	echo "Processing sample $sample_name....."
	echo "Counting $2 > $3 conversions....."
	if [[ ! -d "./$2" ]]; then
		mkdir "./$2"
	fi
	cd $2
	## check if the conversion was processed before
	if [[ -f "counts$2$3.perGene.$sample_name.txt" ]]; then
		echo "$2 to $3 in $sample_name already processed"
		return 0
	fi 
	#first mate, fw T -
	samtools view -f 99 -b -h -L $bed_path $1 | \
	samtools mpileup -A -B -Q 27 -d 1000000 -f $fasta_path - | \
	${cntr_path}/row_mpile_coverage_plus_"$2""$3".pl | \
	awk -v mainbase="$2" -v loopbase="${prs[$2]}" -v base="$3" '($3==loopbase && $5>0){print $1"\t"($2-1)"\t"$2"\t"mainbase""base"\t"$5"\t-";}' > \
	genome"$2""$3"."$sample_name".bed
	#second mate, re T -
	samtools view -f 147 -b -h -L $bed_path $1 | \
	samtools mpileup -A -B -Q 27 -d 1000000 -f $fasta_path - | \
	${cntr_path}/row_mpile_coverage_plus_"$2""$3".pl | \
	awk -v mainbase="$2" -v loopbase="${prs[$2]}" -v base="$3" '($3==loopbase && $5>0){print $1"\t"($2-1)"\t"$2"\t"mainbase""base"\t"$5"\t-";}' >> \
	genome"$2""$3"."$sample_name".bed
	#first mate, re T +
	samtools view -f 83 -b -h -L $bed_path $1 | \
	samtools mpileup -A -B -Q 27 -d 1000000 -f $fasta_path - | \
	${cntr_path}/row_mpile_coverage_plus_"$2""$3".pl | \
	awk -v mainbase="$2" -v loopbase="$2" -v base="$3" '($3==loopbase && $5>0){print $1"\t"($2-1)"\t"$2"\t"mainbase""base"\t"$5"\t+";}' >> \
	genome"$2""$3"."$sample_name".bed
	#second mate, fw +
	samtools view -f 163 -b -h -L $bed_path $1 | \
	samtools mpileup -A -B -Q 27 -d 1000000 -f $fasta_path - | \
	${cntr_path}/row_mpile_coverage_plus_"$2""$3".pl | \
	awk -v mainbase="$2" -v loopbase="$2" -v base="$3" '($3==loopbase && $5>0){print $1"\t"($2-1)"\t"$2"\t"mainbase""base"\t"$5"\t+";}' >> \
	genome"$2""$3"."$sample_name".bed
	
	awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6,$5}' genome"$2""$3"."$sample_name".bed | \
	awk '{a[$1]+=$2;}END{for(i in a)print i"\t"a[i];}' | \
	sed 's/_/\t/g' | \
	awk '{print $1,$2,$3,$4,$7,$6}' OFS="\t" > uniq.genome"$2""$3"."$sample_name".bed

	sort -k1,1 -k2,2n uniq.genome"$2""$3"."$sample_name".bed > sorted.uniq.genome"$2""$3"."$sample_name".bed

	intersectBed -s -wa -wb -sorted -a $bed_path -b sorted.uniq.genome"$2""$3"."$sample_name".bed > "$2""$3".gfeat."$sample_name".bed
	awk '{print $4,$8,$12,$13,$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$8}' OFS="\t" "$2""$3".gfeat."$sample_name".bed | \
	awk '!seen[$5]++' | \
	awk '{a[$1"_"$2]+=$4;}END{for(i in a)print i"\t"a[i];}' > counts"$2""$3".perGene."$sample_name".txt

	# report
	e_time="$(date -u +%s)"
	el_s=$(($e_time - $s_time))
	el_time=$(date --date='@'$el_s +%H:%M:%S)
	echo "Counted $2 > $3 conversions in sample $sample_name"
	echo "Time elapsed: $el_time"
}

# declase a hash of bpairs
declare -A prs=( ["T"]="A" ["C"]="G" ["G"]="C" ["A"]="T" )
# start
# set total start time
ts_time="$(date -u +%s)"
# relocate to the output directory
cd $results_path
# array of bases to loop over
ubloop=( "T" 
	 	 "C"
	 	 "G"
	 	 "A" )
# loop over samples
for x in ${path_to_bams[@]}; do
	# loop over reference bases
	for y in ${ubloop[@]}; do
		# create array of substitutions and drop the reference base from this array
		lbloop=( "T" 
				 "C"
	 	 		 "G"
	 	 		 "A" )
		for j in "${!lbloop[@]}"; do
			if [[ "${lbloop[$j]}" = "$y" ]]; then
				unset 'lbloop[$j]'
			fi
		done
		# loop over subtitutions bases
		for z in ${lbloop[@]}; do
			## for timepoints after 0h only do T -> C
			if echo $x | grep -q -v T0 && [ $y != "T" ] && [ $z != "C" ]; then
				continue
			fi
			(
				countc $x $y $z
			) &
			# allow parallel execution of only N jobs
			if [[ $(jobs -r -p | wc -l) -gt $((N - 1)) ]]; then
        		# wait a batch to finish
        		wait
    		fi
		done
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
