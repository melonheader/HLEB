#!/usr/bin/bash
### INFO: Count genome T -> C conversions from PE Slamseq
### DATE: 11.08.2020
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
cntr_path=$basepath"/scripts/row_mpile_coverage_plus_TC.pl"
wd_path=$basepath"/mapped_reads"
cnv_path=$basepath"/conversions"
out_path=$cnv_path"/genomeTC"
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
echo "Purpose: Count T -> C conversions in .bam files from PE-Slamseq"
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
if [[ ! -f "$cntr_path" ]]; then
	echo ""
	echo "$cntr_path"
	echo "^ .pl counter script is not found. exiting....."
	exit 1
fi


#{CONTROL}
# get the list of BAMs
declare -a file_list=()
while IFS= read -r -d '' file; do
	file_list=("${file_list[@]}" "$file")
done < <(find $wd_path -type f -name '*_Aligned.sortedByCoord.out.bam' -print0)
echo "${#file_list[@]} .bam files found"


#{MAIN}
# define main function $1 - $bam_path
countc() {
	# set process start time
	s_time="$(date -u +%s)"
	# get bam basename and sample_name
	bam_file="$(basename -- $1)"
	sample_name="${bam_file%_S*}"
	echo ""
	echo "Processing sample $sample_name....."
	#first mate, fw T -
	samtools view -f 99 -b -h -L $utr_path $1 | \
	samtools mpileup -A -B -Q 27 -d 1000000 -f $gfa_path - | \
	$cntr_path | \
	awk '($3=="A" && $5>0){print $1"\t"($2-1)"\t"$2"\tTC\t"$5"\t-";}' > genomeTC.$sample_name$bed
	#second mate, re T -
	samtools view -f 147 -b -h -L $utr_path $1 | \
	samtools mpileup -A -B -Q 27 -d 1000000 -f $gfa_path - | \
	$cntr_path | \
	awk '($3=="A" && $5>0){print $1"\t"($2-1)"\t"$2"\tTC\t"$5"\t-";}' >> genomeTC.$sample_name$bed
	#first mate, re T +
	samtools view -f 83 -b -h -L $utr_path $1 | \
	samtools mpileup -A -B -Q 27 -d 1000000 -f $gfa_path - | \
	$cntr_path | \
	awk '($3=="T" && $5>0){print $1"\t"($2-1)"\t"$2"\tTC\t"$5"\t+";}' >> genomeTC.$sample_name$bed
	#second mate, fw +
	samtools view -f 163 -b -h -L $utr_path $1 | \
	samtools mpileup -A -B -Q 27 -d 1000000 -f $gfa_path - | \
	$cntr_path | \
	awk '($3=="T" && $5>0){print $1"\t"($2-1)"\t"$2"\tTC\t"$5"\t+";}' >> genomeTC.$sample_name$bed
	
	awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6,$5}' genomeTC.$sample_name$bed | \
	awk '{a[$1]+=$2;}END{for(i in a)print i"\t"a[i];}' | \
	sed 's/_/\t/g' | \
	awk '{print $1,$2,$3,$4,$7,$6}' OFS="\t" > uniq.genomeTC.$sample_name$bed

	sort -k1,1 -k2,2n uniq.genomeTC.$sample_name$bed > sorted.uniq.genomeTC.$sample_name$bed

	intersectBed -s -wa -wb -sorted -a $utr_path -b sorted.uniq.genomeTC.$sample_name$bed > TC.tpu.$sample_name$bed
	awk '{print $7,$8,$12,$13,$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$7}' OFS="\t" TC.tpu.$sample_name$bed | \
	awk '!seen[$5]++' | \
	awk '{a[$1]+=$4;}END{for(i in a)print i"\t"a[i];}' > genomeTC.perGene.$sample_name.txt

	# report
	e_time="$(date -u +%s)"
	el_s=$(($e_time - $s_time))
	el_time=$(date --date='@'$el_s +%H:%M:%S)
	echo "$sample_name processed!"
	echo "Time elapsed: $el_time"
}

# start
# set total start time
ts_time="$(date -u +%s)"
# relocate to the output directory
cd $out_path
# parallelize over batches of 5
N=5
for i in ${file_list[@]}; do
    (
        countc $i
    ) &
    # allow parallel execution of only N jobs
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait a batch to finish
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
echo "${#file_list[@]} files procced"
echo "Total elapsed time: $tel_time"
echo "Done!"