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
### Auxiliary
if [ -z ${path_to_bams+x} ]; then
	echo "No input bam files provided"
	exit 1
fi
if [ -z ${out_path+x} ]; then
	echo "No output dir specified ayy"
	exit 1
fi
if [ -z ${genome+x} ]; then
	echo "No genome name provided"
	exit 1
fi
if [ -z ${genome_annotation+x} ]; then
	echo "No genome_annotation provided; Using bed of exonic regions"
	genome_annotation="exon"
fi
if [ -z ${seq_style+x} ]; then
	seq_style="SE"
fi
if [[ ! -d "$out_path/$experiment_name" ]]; then
	mkdir -p "$out_path/$experiment_name"
fi
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
bed_path=$(find $annot_dir -type f -iname "*$genome_annotation*.bed")
fasta_path=$annot_dir"/genome_fasta.fa"



#{MAIN}
echo "$(date -u +'%d.%m.%Y')"
echo "Starting $experiment_name!"
echo "Selected genome: $genome with $genome_annotation level annotation"
echo "Output dir: $out_path/$experiment_name"
# set total start time
ts_time="$(date -u +%s)"
## Count bases per genomic feature
bash "$source_dir"/"$seq_style".genomeN.sh \
	-i ${path_to_bams[@]} -o "$out_path/$experiment_name" \
	-b $bed_path -f $fasta_path -n $n_cores > \
	"$out_path"/"$experiment_name"/"$seq_style".genomeN.Log &
wait
## Count base conversions per genomic features
bash "$source_dir"/"$seq_style".genomeNN.sh \
	-i ${path_to_bams[@]} -o "$out_path/$experiment_name" \
	-b $bed_path -f $fasta_path -n $n_cores > \
	"$out_path"/"$experiment_name"/"$seq_style".genomeNN.Log &
wait
## Count SNPs per genomic feature
bash "$source_dir"/genomeSNP.sh \
	-i ${path_to_bams[@]} -o "$out_path/$experiment_name" \
	-b $bed_path -f $fasta_path -v "/local/artem/Scripts/varscan/VarScan.v2.3.9.jar" \
	-n $n_cores > \
	"$out_path"/"$experiment_name"/genomeSNP.Log &
wait
# time report
te_time="$(date -u +%s)"
tel_s=$(($te_time - $ts_time))
tel_time=$(date --date='@'$tel_s +%H:%M:%S)
echo ""
echo "Total elapsed time: $tel_time"
echo "Done!"