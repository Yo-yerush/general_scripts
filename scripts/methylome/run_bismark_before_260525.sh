#!/bin/bash

####################

# download TAIR10 reference genome:
# cd /PATH/TO
# wget -O TAIR10_chr_all.fas.gz https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz

####################

### 'samples_table.txt' example:
### tab delimiter table of sample name and full path to .fastq files
#mut_1 PATH/TO/aaa_R1.fastq
#mut_1 PATH/TO/aaa_R2.fastq
#mut_2	PATH/TO/bbb_R1.fastq
#mut_2 PATH/TO/bbb_R2.fastq
#wt_1  PATH/TO/ccc_R1.fastq
#wt_1  PATH/TO/ccc_R2.fastq
#wt_2  PATH/TO/ddd_R1.fastq
#wt_2  PATH/TO/ddd_R2.fastq

####################

###  Bismark pipeline ###
sample_table="./yo_test.txt"
genome_file_full_path=/PATH/TO/TAIR10_chr_all.fa.gz # default

kepp_just_CX_report=false
sorted_bam_file=true # only if 'kepp_just_CX_report=false'
n_cores=30 # default=10

output_path=./bismark_results # default
output_suffix="wgbs_$(date +%d%m%y)" # default

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

# read sample names and fastq file paths as array
mapfile -t sample_name < <(awk '!seen[$1]++ {print $1}' "$sample_table")
mapfile -t R1_fastq_path < <(awk '$2 ~ /_R1/ {print $2}' "$sample_table")
mapfile -t R2_fastq_path < <(awk '$2 ~ /_R2/ {print $2}' "$sample_table")

paired_end_sequence=true
(( ${#R2_fastq_path[@]} == 0 )) && paired_end_sequence=false
echo "paired-end sequence: $paired_end_sequence"

####################

# n-cores for 'bismark_methylation_extractor'
if [ $n_cores -gt 2 ]; then
      n_cores_2=$((n_cores / 3))
else
      n_cores_2=1
fi

####################

ori_path=$(pwd)
mkdir $output_path

# tmp file for analysis
mkdir $output_path/tmp
cd $output_path/tmp

# index genom
mkdir -p $output_path/genome_indx
cp $genome_file_full_path $output_path/genome_indx
bismark_genome_preparation --verbose $output_path/genome_indx

####################

for ((u=0; u<${#sample_name[@]}; u++)); do
    i="${sample_name[$u]}"
    R1_i="${R1_fastq_path[$u]}"
    R2_i="${R2_fastq_path[$u]:-}"
    
    echo "Processing sample: $i"
    mkdir -p "$output_path/$i"

      ### mapping to genome
      if [[ "$peired_end_sequence" == "false" ]]; then
            echo "mapping to genome for single-end sequence"
            echo "read1 file: $R1_i"
            bismark --bowtie2 -p "$n_cores" $output_path/genome_indx "$R1_i" -o $output_path/"$i" --basename "$i"_$output_suffix
      else
            echo "mapping to genome for peired-end sequence"
            echo "read1 file: $R1_i"
            echo "read2 file: $R2_i"
            bismark --bowtie2 -p "$n_cores" $output_path/genome_indx -1 "$R1_i" -2 "$R2_i" -o $output_path/"$i" --basename "$i"_$output_suffix
      fi

      ### methylation calling
      echo "methylation calling"
      mkdir $output_path/"$i"/methylation_extractor
      bismark_methylation_extractor --CX --parallel "$n_cores_2" --buffer_size 10G --genome_folder $output_path/genome_indx -o $output_path/"$i"/methylation_extractor $output_path/"$i"/"$i"_$output_suffix*.bam

      if [[ "$kepp_just_CX_report" == "true" ]]; then
            echo "keep just 'CX_report' file"
            mv $output_path/"$i"/methylation_extractor/*.CX_report.txt $output_path
            rm -r $output_path/"$i"

      else
            echo "keep output '.bam' and '.cov' files"
            if [[ "$sorted_bam_file" == "true" ]]; then
                  echo "sort bam files (can use in IGV softwar)"
                  samtools sort $output_path/"$i"/"$i"_$output_suffix*.bam -o $output_path/"$i"/"$i"_"$output_suffix"_sorted.bam
                  samtools index $output_path/"$i"/"$i"_"$output_suffix"_sorted.bam
            fi
      fi

      ((u++))
      echo "Completed sample: $i"
      echo "-----------------------------------"
      echo ""
done

if [[ "$kepp_just_CX_report" == "true" ]]; then
      rm -r $output_path/genome_indx
fi

cd $ori_path
rm -r $output_path/tmp
