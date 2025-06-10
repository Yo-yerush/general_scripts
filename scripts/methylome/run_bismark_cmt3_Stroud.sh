#!/bin/bash

###  Bismark pipeline ###
samples_name=("SRR534177" "SRR534193" "SRR534209")
fastq_file_suffix="fastq.gz"
R1_suffix="_1"
R2_suffix="_2"
output_path=/home/yoyerush/yo/Methylome.At_project/paper_files/Stroud_et_al_data/bismark_results
fastq_directory=/home/yoyerush/yo/Methylome.At_project/paper_files/Stroud_et_al_data/fastq_files
output_suffix=cmt3_Stroud
genome_file_full_path=/home/yoyerush/yo/Methylome.At_project/paper_files/TAIR10_chr_all.fa.gz
n_cores=30

peired_end_sequence=false

kepp_just_CX_report=true
sorted_bam_file=false # only if 'kepp_just_CS_report=false'

# n-cores for 'bismark_methylation_extractor'
if [ $n_cores -gt 2 ]; then
      n_cores_2=$((n_cores / 3))
else
      n_cores_2=1
fi

####################
# download TAIR10 reference genome at:
# https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
####################

#### add sample_vec instead of numbers in the for loop
mkdir $output_path

# tmp file for analysis
ori_path=pwd
mkdir $output_path/tmp
cd $output_path/tmp

# index genom
mkdir $output_path/genome_indx
cp $genome_file_full_path $output_path/genome_indx
bismark_genome_preparation --verbose $output_path/genome_indx

for i in "${samples_name[@]}"
do
      # mapping to genome for single-end sequence
      mkdir $output_path/"$i"
      if [[ "$peired_end_sequence" == "false" ]]; then
            bismark --bowtie2 -p "$n_cores" $output_path/genome_indx $fastq_directory/"$i"."$fastq_file_suffix" -o $output_path/"$i" --basename "$i"_$output_suffix
      else
            bismark --bowtie2 -p "$n_cores" $output_path/genome_indx -1 $fastq_directory/"$i"*"$R1_suffix"."$fastq_file_suffix" -2 $fastq_directory/"$i"*"$R2_suffix"."$fastq_file_suffix" -o $output_path/"$i" --basename "$i"_$output_suffix
      fi
      
      # methylation calling
      if [[ "$kepp_just_CX_report" == "true" ]]; then
            # keep just 'CX_report' file
            bismark_methylation_extractor --CX --parallel "$n_cores_2" --buffer_size 10G --genome_folder $output_path/genome_indx -o $output_path/"$i"/methylation_extractor $output_path/"$i"/"$i"_$output_suffix*.bam

            mv $output_path/"$i"/methylation_extractor/*.CX_report.txt $output_path
            rm -r $output_path/"$i"
            rm -r $output_path/genome_indx
      else
            # keep output '.bam' and '.cov' files
            mkdir $output_path/"$i"/methylation_extractor
            bismark_methylation_extractor --parallel "$n_cores_2" --buffer_size 10G --genome_folder $output_path/genome_indx -o $output_path/"$i"/methylation_extractor $output_path/"$i"/"$i"_$output_suffix*.bam

            # sort bam files (can use in IGV softwar)
            if [[ "$sorted_bam_file" == "true" ]]; then
                  samtools sort $output_path/"$i"/"$i"_$output_suffix*.bam -o $output_path/"$i"/"$i"_"$output_suffix"_sorted.bam
                  samtools index $output_path/"$i"/"$i"_"$output_suffix"_sorted.bam
            fi
      fi
done

cd $ori_path
rm -r $output_path/tmp