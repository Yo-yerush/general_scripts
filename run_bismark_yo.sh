#!/bin/bash
###############################################################################
# YO - 260525
# Bismark WGBS pipeline
#
# example:
# download TAIR10 reference genome:
# cd /PATH/TO
# wget -O TAIR10_chr_all.fas.gz https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
#
# Run:
#   ./run_bismark_yo.sh -s samples_table.txt -g TAIR10_chr_all.fa.gz
#
# Usage
# -----
#   run_bismark_yo.sh  -s|--samples  FILE          (required)
#                   -g|--genome   FASTA.gz      (required)
#                   -o|--outdir   DIR           [default: ./bismark_results]
#                   -p|--suffix   STRING        [default: wgbs_<date (d-m-y)>]
#                   -t|--threads  N             [default: 10]
#                   --keep-cx                   (keep only *.CX_report.txt)
#                   --sort-bam                  (produce sorted & indexed BAM)
#                   -h|--help
#
# Options
# -------
#   -s, --samples   Tab-delimited two-column file: sample-name <TAB> fastq-path
#   -g, --genome    Compressed FASTA of the reference genome (will be indexed)
#   -o, --outdir    Where to write all results and temporary files
#   -p, --suffix    Basename suffix applied to all Bismark outputs
#   -t, --threads   Number of CPU cores to use (Bowtie2; extractor uses ⌊N/3⌋ )
#   --keep-cx       Delete everything except *.CX_report.txt (saves space)
#   --sort-bam      After extraction, sort & index BAM for IGV viewing
#   -h, --help      Show this message and exit
###############################################################################

### default values
sample_table=
genome_file_full_path=
output_path="./bismark_results"
output_suffix="wgbs_$(date +%d%m%y)"
n_cores=10
keep_cx=false
sort_bam=false

while [[ $# -gt 0 ]]; do
      case $1 in
      -s | --samples)
            sample_table=$2
            shift 2
            ;;
      -g | --genome)
            genome_file_full_path=$2
            shift 2
            ;;
      -o | --outdir)
            output_path=$2
            shift 2
            ;;
      -p | --prefix)
            output_suffix=$2
            shift 2
            ;;
      -t | --threads)
            n_cores=$2
            shift 2
            ;;
      --keep-cx)
            keep_cx=true
            shift
            ;;
      --sort-bam)
            sort_bam=true
            shift
            ;;
      -h | --help)
            sed -n '1,/^###############################################################################$/p' "$0"
            exit 0
            ;;
      *)
            echo "Unknown option: $1"
            exit 1
            ;;
      esac
done

####################

# check if 'sample_table' file exists
if [[ ! -f "$sample_table" ]]; then
      echo "Error: Sample table file '$sample_table' does not exist."
      exit 1
fi
# check if 'genome_file_full_path' file exists
if [[ ! -f "$genome_file_full_path" ]]; then
      echo "Error: Genome file '$genome_file_full_path' does not exist."
      exit 1
fi

### read sample names and fastq file paths as array
mapfile -t sample_name < <(awk '!seen[$1]++ {print $1}' "$sample_table")
mapfile -t R1_fastq_path < <(awk '$2 ~ /(_R1|_1\.f)/ {print $2}' "$sample_table")
mapfile -t R2_fastq_path < <(awk '$2 ~ /(_R2|_2\.f)/ {print $2}' "$sample_table")

paired_end_sequence=true
((${#R2_fastq_path[@]} == 0)) && paired_end_sequence=false

####################

# n-cores for 'bismark_methylation_extractor'
#if [ $n_cores -gt 2 ]; then
#      n_cores_2=$((n_cores / 3))
#else
#      n_cores_2=1
#fi

####################

ori_path=$(pwd)
mkdir $output_path
cd $output_path
output_path=$(pwd)

### tmp file for analysis
mkdir $output_path/tmp
cd $output_path/tmp

### Generate log file with a timestamp
log_file="../${output_suffix}_.log"
echo "**  $(date +"%d-%m-%y %H:%M")" >"$log_file"
echo "**  samples: ${sample_name[@]}" >>"$log_file"
#echo "paired-end sequence: $paired_end_sequence" >>"$log_file"

echo "" >>"$log_file"

### index genom
genome_b_name=$(basename "$genome_file_full_path")
genome_new_path=$output_path/genome_indx/$genome_b_name
echo "indexing genome file: '$genome_b_name'" >> "$log_file"
mkdir -p $output_path/genome_indx
cp $genome_file_full_path $output_path/genome_indx
# Rename genome file if necessary
if [[ "$genome_b_name" == *.fas ]]; then
      mv "$genome_new_path" "${genome_new_path%.fas}.fa"
      echo "* rename genome file: '$(basename $genome_file_full_path)'" >> "$log_file"
elif [[ "$genome_new_path" == *.fas.gz ]]; then
      mv "$genome_new_path" "${genome_new_path%.fas.gz}.fa.gz"
      echo "* rename genome file: '$(basename $genome_file_full_path)'" >> "$log_file"
fi
bismark_genome_preparation $output_path/genome_indx

echo "" >>"$log_file"

####################

for ((u = 0; u < ${#sample_name[@]}; u++)); do
      i="${sample_name[$u]}"
      R1_i="${R1_fastq_path[$u]}"
      R2_i="${R2_fastq_path[$u]:-}"

      echo "Processing sample: $i" >>"$log_file"
      mkdir -p "$output_path/$i"

      ### mapping to genome
      if [[ "$peired_end_sequence" == "false" ]]; then
            Rs_type="se"
            echo "mapping to genome for single-end sequence" >>"$log_file"
            echo "read1 file: $R1_i" >>"$log_file"
            bismark --bowtie2 --parallel "$n_cores" $output_path/genome_indx "$R1_i" -o $output_path/"$i" # --basename "$i"_$output_suffix
      else
            Rs_type="pe"
            echo "mapping to genome for peired-end sequence:" >>"$log_file"
            echo "* read1 file: $(basename "$R1_i")" >>"$log_file"
            echo "* read2 file: $(basename "$R2_i")" >>"$log_file"
            bismark --bowtie2 --parallel "$n_cores" $output_path/genome_indx -1 "$R1_i" -2 "$R2_i" -o $output_path/"$i" # --basename "$i"_$output_suffix
      fi

      ### methylation calling
      echo "" >>"$log_file"
      echo "methylation calling" >>"$log_file"
      mkdir $output_path/"$i"/methylation_extractor
      bismark_methylation_extractor --CX --parallel "$n_cores" --genome_folder $output_path/genome_indx -o $output_path/"$i"/methylation_extractor $output_path/"$i"/*_bt2_"$Rs_type".bam #"$i"_$output_suffix*.bam

      if [[ "$keep_cx" == "true" ]]; then
            echo "* keep just 'CX_report' file" >>"$log_file"
            mv $output_path/"$i"/methylation_extractor/*.CX_report.txt $output_path/"$i"/methylation_extractor/"$i".CX_report.txt # change name
            mv $output_path/"$i"/methylation_extractor/"$i".CX_report.txt $output_path
            rm -r $output_path/"$i"

      else
            echo "* keep output '.bam' and '.cov' files" >>"$log_file"
            if [[ "$sort_bam" == "true" ]]; then
                  echo "* sort bam files (can use in IGV softwar)" >>"$log_file"
                  samtools sort $output_path/"$i"/"$i"_$output_suffix*.bam -o $output_path/"$i"/"$i"_"$output_suffix"_sorted.bam
                  samtools index $output_path/"$i"/"$i"_"$output_suffix"_sorted.bam
            fi
      fi

      ((u++))
      echo "Completed sample: $i" >>"$log_file"
      echo "" >>"$log_file"
      echo "-----------------------------------" >>"$log_file"
      echo "" >>"$log_file"
done

if [[ "$keep_cx" == "true" ]]; then
      rm -r $output_path/genome_indx
fi

echo "**  $(date +"%d-%m-%y %H:%M")" >>"$log_file"
cd $ori_path
rm -r $output_path/tmp