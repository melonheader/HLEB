#!/usr/bin/bash
### INFO: Collapse GTFs into BED files for a genomic feature
### DATE: 02.08.2022
### AUTHOR: Artem Baranovskii

gtf=$1
feature=$2
out_path=$3
out_file=$(basename -- $gtf .gtf)_"$feature"_collapsed.bed

if [[ $feature == "exon" ]]; then
    offset=2
else 
    offset=0
fi

cd $out_path
cat $gtf | \
awk -v genomic_feature=$feature -v c_offset=$offset 'BEGIN{OFS="\t";} $3==genomic_feature {print $1,$4-1,$5,$10,"1",$7,genomic_feature,$(18 + c_offset)}' | \
tr -d '"|;' | \
sortBed | \
groupBy -g 4,6,7,8 -c 1,2,3 -o distinct,min,max | \
awk -v OFS='\t' '{print $5, $6, $7, $1, 1, $2, $3, $4}' > "$out_path"/"$out_file"
#mergeBed -s -c 4,5,6,7,8 -o distinct,distinct,distinct,distinct,distinct > "$out_path"/"$out_file" 
# ^ resulted in merging of overlapping regions that are laborous to disentangle later, use groupBy instead to keep them separate