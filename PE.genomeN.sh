#!/usr/bin/bash
### INFO: Count genome N for conversion rate estimation from PE Slamseq
### DATE: 04.09.2020
### AUTHOR: Miha Milek, Artem Baranovskii


#{SETUP}
set -e
# set variables
# ------------------------------ #
## Global
project="Stability"
exp="Slamseq_caf1DN"
genome="mm10"
# ----------------------------------------------------------------- #
## Paths
basepath="/local/artem/Projects/"$project"/Data/SequenceData/"$exp
utr_path="/local/shared_scripts/database/"$genome"/sorted.Mus_musculus.GRCm38.96.3utr.bed"
gfa_path="/local/shared_scripts/database/"$genome"/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
cntr_path="/local/artem/Scripts/counters"
wd_path=$basepath"/mapped_reads"
cnv_path=$basepath"/conversions"
out_path=$cnv_path"/genomeN"
### Extensions
bam=".bam"
sam=".sam"
prefix="sorted.allReads."
suffix="_Aligned.sortedByCoord.out"
bed=".bed"
txt=".txt"
pileup=".pileup"


# ----------------------------------------------------------------- #
# Echo settings
echo "Settings:"
echo "Project: $project"
echo "Experiment: $exp"
echo "Purpose: Count genome N for conversion rate estimation in .bam files from PE-Slamseq"
echo ""
# check for fasta files
if [[ ! -d "$wd_path" ]]; then
	echo ""
	echo "$wd_path"
	echo "^ folder with .bam files is not found. exiting....."
	exit 1
fi
# check for out dir
if [[ ! -d "$cnv_path" ]]; then
	mkdir "$cnv_path"
fi
if [[ ! -d "$out_path" ]]; then
	mkdir "$out_path"
fi
# check for .pl counter
if [[ ! -d "$cntr_path" ]]; then
	echo ""
	echo "$cntr_path"
	echo "^ .pl counters folder is not found. exiting....."
	exit 1
fi


#{CONTROL}
# get the list of BAMs
declare -a file_list=()
while IFS= read -r -d '' file; do
	file_list=("${file_list[@]}" "$file")
done < <(find $wd_path -type f -name '*_Aligned.sortedByCoord.out.bam' -print0)
echo "${#file_list[@]} .bam files found"
# keep only timepoint 0
file_sel=$( for i in ${file_list[@]}; do echo $i; done | grep T0 )
echo ""
echo "${#file_sel[@]} .bam files selected"
# declase a hash of bpairs
declare -A bloop=( ["T"]="t" ["C"]="c" ["G"]="g" ["A"]="a" )


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
	#first mate, fw T -
	samtools view -f 99 -b -h -L $utr_path $1 | \
	samtools mpileup -A -B -Q 1 -d 1000000 -f $gfa_path - | \
	awk -v base="$2" -v lbase="${bloop[$2]}" '($3==lbase || $3==base){print $1,$2,$3,$4,$5,"-"}' OFS="\t" > \
	cov$2.$sample_name$pileup
	#second mate, re T -
	samtools view -f 147 -b -h -L $utr_path $1 | \
	samtools mpileup -A -B -Q 1 -d 1000000 -f $gfa_path - | \
	awk -v base="$2" -v lbase="${bloop[$2]}" '($3==lbase || $3==base){print $1,$2,$3,$4,$5,"-"}' OFS="\t" >> \
	cov$2.$sample_name$pileup
	#first mate, re T +
	samtools view -f 83 -b -h -L $utr_path $1 | \
	samtools mpileup -A -B -Q 1 -d 1000000 -f $gfa_path - | \
	awk -v base="$2" -v lbase="${bloop[$2]}" '($3==lbase || $3==base){print $1,$2,$3,$4,$5,"+"}' OFS="\t" >> \
	cov$2.$sample_name$pileup
	#second mate, fw +
	samtools view -f 163 -b -h -L $utr_path $1 | \
	samtools mpileup -A -B -Q 1 -d 1000000 -f $gfa_path - | \
	awk -v base="$2" -v lbase="${bloop[$2]}" '($3==lbase || $3==base){print $1,$2,$3,$4,$5,"+"}' OFS="\t" >> \
	cov$2.$sample_name$pileup
	
	awk '{gsub(">","<",$5); print}' cov$2.$sample_name$pileup | \
	awk '{print $1"_"$2"_"$3"_"$6,$4-gsub(/</,"",$5)}' OFS="\t" | \
	awk '($2>0)' | \
	awk '{a[$1]+=$2;}END{for(i in a)print i"\t"a[i];}' | \
	sed 's/_/\t/g' | \
	awk '{print $1,$2,$2+1,$3,$5,$4}' OFS="\t" > cov$2.$sample_name$bed

	sort -k1,1 -k2,2n cov$2.$sample_name$bed > sorted.cov$2.$sample_name$bed

	intersectBed -s -wa -wb -sorted -a $utr_path -b sorted.cov$2.$sample_name$bed > cov$2.tpu.$sample_name$bed
	awk '{print $7,$8,$12,$13,$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$7}' OFS="\t" cov$2.tpu.$sample_name$bed | \
	awk '!seen[$5]++' |  \
	awk '{a[$1]+=$4;}END{for(i in a)print i"\t"a[i];}' > counts$2.perGene.$sample_name$txt

	# report
	e_time="$(date -u +%s)"
	el_s=$(($e_time - $s_time))
	el_time=$(date --date='@'$el_s +%H:%M:%S)
	echo "genome background $2 in sample $sample_name processed"
	echo "Time elapsed: $el_time"
}


# start
# set total start time
ts_time="$(date -u +%s)"
# relocate to the output directory
cd $out_path
# parallelize over batches of 6
N=6
# loop over samples
for x in ${file_sel[@]}; do
	# loop over background bases
	for y in ${!bloop[@]}; do
		(
			countb $x $y
		) &
		# allow parallel execution of N jobs
		if [[ $(jobs -r -p | wc -l) -gt $((N - 1))  ]]; then
			# wait for a batch to finish
			wait
		fi
	done
done

