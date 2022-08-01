#!/usr/bin/bash
### INFO: Count genome Ts from SE Slamseq
### DATE: 11.08.2020
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
wd_path=$basepath"/mapped_reads"
cnv_path=$basepath"/conversions"
out_path=$cnv_path"/genomeT"
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
echo "Purpose: Count genome Ts in the provided .bam files from SE-Slamseq"
# check for fasta files
if [[ ! -d "$wd_path" ]]; then
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


#{CONTROL}
# get the list of BAMs
declare -a file_list=()
while IFS= read -r -d '' file; do
	file_list=("${file_list[@]}" "$file")
done < <(find $wd_path -type f -name '*_Aligned.sortedByCoord.out.bam' -print0)
echo "${#file_list[@]} .bam files found"


#{MAIN}
# total elapsed time
ts_time="$(date -u +%s)"
# relocate to the output directory
cd $out_path
# iterate
for bam_path in "${file_list[@]}"
do
	s_time="$(date -u +%s)"
	# get bam basename and sample_name
	bam_file="$(basename -- $bam_path)"
	sample_name="${bam_file%_S*}"
	echo "Processing sample $sample_name....."
	# +
	samtools view -f 16 -b -h -L $utr_path $bam_path | \
	samtools mpileup -A -B -Q 1 -d 1000000 -f $gfa_path - | \
	awk '($3=="t")||($3=="T"){print $1,$2,$3,$4,$5,"-"}' OFS="\t" > covT.$sample_name$pileup
	# -
	samtools view -F 16 -b -h -L $utr_path $bam_path | \
	samtools mpileup -A -B -Q 1 -d 1000000 -f $gfa_path - | \
	awk '($3=="t")||($3=="T"){print $1,$2,$3,$4,$5,"-"}' OFS="\t" >> covT.$sample_name$pileup

	
	awk '{gsub(">","<",$5); print}' covT.$sample_name$pileup | \
	awk '{print $1"_"$2"_"$3"_"$6,$4-gsub(/</,"",$5)}' OFS="\t" | \
	awk '($2>0)' | \
	awk '{a[$1]+=$2;}END{for(i in a)print i"\t"a[i];}' | \
	sed 's/_/\t/g' | \
	awk '{print $1,$2,$2+1,$3,$5,$4}' OFS="\t" > covT.$sample_name$bed

	sort -k1,1 -k2,2n covT.$sample_name$bed > sorted.covT.$sample_name$bed

	intersectBed -s -wa -wb -sorted -a $utr_path -b sorted.covT.$sample_name$bed > covT.tpu.$sample_name$bed
	awk '{print $7,$8,$12,$13,$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$7}' OFS="\t" covT.tpu.$sample_name$bed | \
	awk '!seen[$5]++' |  \
	awk '{a[$1]+=$4;}END{for(i in a)print i"\t"a[i];}' > countsT.perGene.$sample_name$txt

	# report
	e_time="$(date -u +%s)"
	el_s=$(($e_time - $s_time))
	el_time=$(date --date='@'$el_s +%H:%M:%S)
	echo "$sample_name processed!"
	echo "Time elapsed: $el_time"
done

# t report
te_time="$(date -u +%s)"
tel_s=$(($te_time - $ts_time))
tel_time=$(date --date='@'$tel_s +%H:%M:%S)
echo "${#file_list[@]} files procced"
echo "Total elapsed time: $tel_time"
echo "Done!"