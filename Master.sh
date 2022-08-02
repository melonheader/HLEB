#!/usr/bin/env bash
### INFO: Process Slamseq data
### DATE: 01.08.2022
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
			echo "-g, --genome			  name of the genome <mm10|GRCh38|rnor6>"
			echo "-a, --genome_annotation name of the genome annotation to use <three_prime_utr|exon>"
			echo "-e, --experiment_name   name of the current run, DD/MM/YY if not specified"
            exit 0
        ;;
        -i|--path_to_bams)
            shift 
            if test $# -gt 0; then
                path_to_bams=$1
            else
                echo "No input bam files provided"
				exit 1
            fi
            shift
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
		-g|--genome)
			shift
			if test $# -gt 0; then
				genome=$1
			else 
				echo "No genome name provided"
				exit 1
			fi
			shift
		;;
		-a|--genome_annotation)
			shift 
			if test $# -gt 0; then
				genome_annotation=$1
			else
				echo "No genome_annotation provided\nUsing bed of exonic regions"
				genome_annotation="exon"
			fi
			shift
		;;
		-paired)
			seq_style="PE"
        ;;
		-e|--experiment_name)
			shift 
			if test $# -gt 0; then
				experiment_name=$1
			else 
				experiment_name="$(date -u +'%d.%m.%Y')_Run"
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
### Annotation 
annot_dir="/local/shared_scripts/database/"$genome"/slam_annot"
bed_path=$(find $annot_dir -type f -iname "*'$genome_annotation'_collapsed.bed")
fasta_path=$annot_dir"/genome_fasta.fa"
### Auxiliary
if [ -z ${seq_style+x} ]; then
	seq_style="SE"
fi


#{MAIN}
## Count bases per genomic feature
bash "$source_dir"/"$seq_style".genomeN.sh \
	-i $path_to_bams -o $out_path \
	-b $bed_path -f $fasta_path -e $experiment_name -n 6 > \
	"$out_path"/"$seq_style".genomeN.Log
## Count base conversions per genomic features
bash "$source_dir"/"$seq_style".genomeNN.sh \
	-i $path_to_bams -o $out_path \
	-b $bed_path -f $fasta_path -e $experiment_name -n 6 > \
	"$out_path"/"$seq_style".genomeNN.Log
## Count SNPs per genomic feature
bash "$source_dir"/genomeSNP.sh \
	-i $path_to_bams -o $out_path \
	-b $bed_path -f $fasta_path -v "/local/artem/Scripts/varscan/VarScan.v2.3.9.jar" \
	-e $experiment_name -n 6 > \
	"$out_path"/genomeSNP.Log