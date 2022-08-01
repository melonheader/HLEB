#!/usr/bin/bash
### INFO: Count genome N -> N conversions from SE Slamseq
### DATE: 03.09.2020
### AUTHOR: Miha Milek, Artem Baranovskii


#{SETUP}
set -e
# set variables
# ------------------------------ #
## Global
project="Stability"
exp="PCN.Slamseq_dnCaf1"
genome="mm10"
# ----------------------------------------------------------------- #
## Paths
basepath="/local/artem/Projects/"$project"/Data/SequenceData/"$exp
utr_path="/local/shared_scripts/database/"$genome"/sorted.Mus_musculus.GRCm38.96.3utr.bed"
gfa_path="/local/shared_scripts/database/"$genome"/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
cntr_path="/local/artem/Scripts/counters"
wd_path=$basepath"/mapped_reads"
cnv_path=$basepath"/conversions"
out_path=$cnv_path"/genomeNN"
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
echo "Purpose: Count N -> N conversions in .bam files from PE-Slamseq"
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
file_sel=$( for i in ${file_list[@]} ; do echo $i ; done | grep "0hrs" )
echo ""
echo "${#file_sel[@]} .bam files selected"
# declase a hash of bpairs
declare -A prs=( ["T"]="A" ["C"]="G" ["G"]="C" ["A"]="T" )


#{MAIN}
# define main function $1 - bam_path, $2 - main base, $3 - rev base
countc() {
	# set process start time
	s_time="$(date -u +%s)"
	# get bam basename and sample_name
	bam_file="$(basename -- $1)"
	sample_name="${bam_file%_S*}"
	echo ""
	echo "Processing sample $sample_name....."
	echo "Counting $2 > $3 conversions....."
	if [[ ! -d "./$2" ]]; then
		mkdir "./$2"
	fi
	cd $2

	# +
	samtools view -f 16 -b -h -L $utr_path $1 | \
	samtools mpileup -A -B -Q 27 -d 1000000 -f $gfa_path - | \
	${cntr_path}/row_mpile_coverage_plus_$2$3.pl | \
	awk -v mainbase="$2" -v loopbase="${prs[$2]}" -v base="$3" '($3==loopbase && $5>0){print $1"\t"($2-1)"\t"$2"\t"mainbase""base"\t"$5"\t-";}' > \
	genome$2$3.$sample_name$bed

	# -
	samtools view -F 16 -b -h -L $utr_path $1 | \
	samtools mpileup -A -B -Q 27 -d 1000000 -f $gfa_path - | \
	${cntr_path}/row_mpile_coverage_plus_$2$3.pl | \
	awk -v mainbase="$2" -v loopbase="${prs[$2]}" -v base="$3" '($3==loopbase && $5>0){print $1"\t"($2-1)"\t"$2"\t"mainbase""base"\t"$5"\t-";}' >> \
	genome$2$3.$sample_name$bed
	
	# filter unique
	awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6,$5}' genome$2$3.$sample_name$bed | \
	awk '{a[$1]+=$2;}END{for(i in a)print i"\t"a[i];}' | \
	sed 's/_/\t/g' | \
	awk '{print $1,$2,$3,$4,$7,$6}' OFS="\t" > uniq.genome$2$3.$sample_name$bed

	# sort
	sort -k1,1 -k2,2n uniq.genome$2$3.$sample_name$bed > sorted.uniq.genome$2$3.$sample_name$bed

	# quantify conversions per gene
	intersectBed -s -wa -wb -sorted -a $utr_path -b sorted.uniq.genome$2$3.$sample_name$bed > $2$3.tpu.$sample_name$bed
	awk '{print $7,$8,$12,$13,$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$7}' OFS="\t" $2$3.tpu.$sample_name$bed | \
	awk '!seen[$5]++' | \
	awk '{a[$1]+=$4;}END{for(i in a)print i"\t"a[i];}' > genome$2$3.perGene.$sample_name.txt

	# report
	e_time="$(date -u +%s)"
	el_s=$(($e_time - $s_time))
	el_time=$(date --date='@'$el_s +%H:%M:%S)
	echo "Counted $2 > $3 conversions in sample $sample_name"
	echo "Time elapsed: $el_time"
}


# start
# set total start time
ts_time="$(date -u +%s)"
# relocate to the output directory
cd $out_path
# array of bases to loop over
ubloop=( "T" 
	 	 "C"
	 	 "G"
	 	 "A" )
# parallelize over batches of N
N=4
# loop over samples
for x in ${file_sel[@]}; do
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
echo "${#file_sel[@]} files procced"
echo "Total elapsed time: $tel_time"
echo "Done!"
