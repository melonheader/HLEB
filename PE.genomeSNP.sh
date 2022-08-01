#!/bin/bash
#$ -V
#$ -l h_vmem=60G,h_rt=48:0:0
#$ -o /scratch/AG_Landthaler/miha/dyn/total/star_endtoend/genomeSNP/
#$ -e /scratch/AG_Landthaler/miha/dyn/total/star_endtoend/genomeSNP/


basepath="/scratch/AG_Landthaler/miha/dyn/total/star_endtoend/sortedAll/"
bam=".bam"
sam=".sam"
prefix="sorted.allReads."
suffix=""
outpath="/scratch/AG_Landthaler/miha/dyn/total/star_endtoend/genomeSNP/"
genome="/scratch/AG_Landthaler/genomes/mm10/Indices/genome.fa"
rscript="/scratch/AG_Landthaler/miha/dyn/TClist.Rscript"
bed=".bed"
minVarFreq="0.8"
minCov="10"


array=( "t0min1" "t15min1" "t20min1" "t30min1" "t45min1" "t60min1" "t0min2" "t15min2" "t20min2" "t30min2" "t45min2" "t60min2" "t60minIAA1" "t60minIAA2" )

cd $outpath

for filename in "${array[@]}"
do


	samtools mpileup -B -A -q 255 -f $genome $basepath$prefix$filename$suffix$bam | java -jar /home/mmilek/varscan/VarScan.v2.3.9.jar  mpileup2snp --strand-filter 0 --output-vcf --min-var-freq $minVarFreq --min-coverage $minCov --variants 1 | awk '($4=="T"&&$5=="C")' | awk '{print $1,$2,$2+1,$4$5,$10}' OFS="\t" | sed 's/:/\t/g' | awk '{print $1,$2,$3,$4,$10,"+"}' OFS="\t" > counts.$filename$bed

	samtools mpileup -B -A -q 255 -f $genome $basepath$prefix$filename$suffix$bam | java -jar /home/mmilek/varscan/VarScan.v2.3.9.jar  mpileup2snp --strand-filter 0 --output-vcf --min-var-freq $minVarFreq --min-coverage $minCov --variants 1 | awk '($4=="A"&&$5=="G")' | awk '{print $1,$2,$2+1,$4$5,$10}' OFS="\t" | sed 's/:/\t/g' | awk '{print $1,$2,$3,$4,$10,"-"}' OFS="\t" >> counts.$filename$bed

	sort -k1,1 -k2,2n counts.$filename$bed > sorted.counts.$filename$bed
	
	intersectBed -s -wa -wb -sorted -a /scratch/AG_Landthaler/miha/star/mm10/exons.gencode.vM14.bed -b sorted.counts.$filename$bed | awk '($6=="+"&&$12=="TC")||($6=="-"&&$12=="AG")'> counts.exons.$filename$bed
	
	awk '{print $7,$8,$12,$13,$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$7}' OFS="\t" counts.exons.$filename$bed | awk '!seen[$5]++' |  awk '{a[$1]+=$4;}END{for(i in a)print i"\t"a[i];}' > countsTC.perGene."$filename".txt

	rm -f counts.$filename$bed
	
done 





