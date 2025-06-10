#!/bin/bash
echo "
               ,@@@@@@@,
       ,,,.   ,@@@@@@/@@,  .oo8888o.
    ,&%%&%&&%,@@@@@/@@@@@@,8888\88/8o
   ,%&\%&&%&&%,@@@\@@@/@@@88\88888/88'
   %&&%&%&/%&&%@@\@@/ /@@@88888\88888'
   %&&%/ %&%%&&@@\ V /@@' '88\8 '/88'
   '&%\ ' /%&'    |.|        \ '|8'
       |o|        | |         | |
       |.|        | |         | |
jgs \\/ ._\//_/__/  ,\_//__\\/.  \_//__/_

"

# configuration
current_date=$(date +%d%m%y)
genome_file_full_path=/home/yoyerush/yo/pomegranate/GCF_007655135.1_ASM765513v2_rna.fna
ann_file_full_path=/home/yoyerush/yo/pomegranate/gene2transcripts_yo.txt
fastq_directory=/home/yaelh/rnaseq_2024_b/Order_582136
output_path=/home/yaelh/rnaseq_2024_b/results_$current_date
output_suffix=pome
plant=pome

to_do_trimming=true
to_do_fastqc=true
is_paired=true

if [ "$to_do_trimming" = true ]; then
new_fastq_directory=$output_path/trimmed_files
else
new_fastq_directory=$fastq_directory
fi

# fastq files names
fastq_file_prefix="S" # before the sample number
fastq_file_suffix_R1="_R1.fastq.gz" # after the sample number
fastq_file_suffix_R2="_R2.fastq.gz" # after the sample number
n_samples=$(seq 1 2) # vector of sample numbers

# .log file
mkdir $output_path
date > $output_path/pipeline_$current_date.log
echo "
configuration:
genome_file_full_path = $genome_file_full_path
ann_file_full_path = $ann_file_full_path
fastq_directory = $fastq_directory
output_path = $output_path
output_suffix = $output_suffix
plant = $plant
to_do_trimming = $to_do_trimming
to_do_fastqc = $to_do_fastqc
is_paired = $is_paired
n_samples = $(echo $n_samples | tr '\n' ' ')

" | tee -a $output_path/pipeline_$current_date.log # print to the console and the .log file

# conda env for trimming and fastqc
source activate fastqc_and_trimming

### trimming
if [ "$to_do_trimming" = true ]; then
    mkdir $output_path/trimmed_files
    cd $output_path/trimmed_files
    for i in $n_samples
    do
    if [ "$is_paired" = true ]; then
        java -jar /usr/bin/trimmomatic.jar PE -threads 20 -phred33 $fastq_directory/${fastq_file_prefix}"$i"${fastq_file_suffix_R1} $fastq_directory/${fastq_file_prefix}"$i"${fastq_file_suffix_R2} S"$i"_R1_paired.fastq.gz S"$i"_R1_unpaired.fastq.gz S"$i"_R2_paired.fastq.gz S"$i"_R2_unpaired.fastq.gz ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10
    else
        java -jar /usr/bin/trimmomatic.jar SE -threads 20 -phred33 $fastq_directory/${fastq_file_prefix}"$i"${fastq_file_suffix_R1} S"$i"_R1.fastq.gz ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10
    fi
    done
    rm $output_path/trimmed_files/*_unpaired.fastq.gz
    echo "trimming: done" >> $output_path/pipeline_$current_date.log
    cd $output_path
fi

### fastqc
if [ "$to_do_trimming" = true ]; then
    mkdir $output_path/qc_files
    cd $output_path/qc_files
    for i in $n_samples
    do
        if [ "$is_paired" = true ]; then
        fastqc -o $output_path/qc_files $new_fastq_directory/S"$i"_R1_paired.fastq.gz
        fastqc -o $output_path/qc_files $new_fastq_directory/S"$i"_R2_paired.fastq.gz
        else
        fastqc -o $output_path/qc_files $new_fastq_directory/S"$i"_R1.fastq.gz
        fi
    done
    rm $output_path/qc_files/*_fastqc.zip
    echo "fastqc: done" >> $output_path/pipeline_$current_date.log
    cd $output_path
fi

echo "" >> $output_path/pipeline_$current_date.log
conda deactivate

### prepare reference
mkdir $output_path/"$plant"_index
rsem-prepare-reference --bowtie2 --transcript-to-gene-map $ann_file_full_path $genome_file_full_path $output_path/"$plant"_index/"$plant"_index

echo "prepare reference: done" >> $output_path/pipeline_$current_date.log

mkdir $output_path/calc_expression
mkdir $output_path/sorted_bam
for i in $n_samples
do
        ### calculate expression
        mkdir $output_path/calc_expression/$output_suffix"$i"
        if [ "$is_paired" = true ]; then
        rsem-calculate-expression -p 20 --bowtie2 --paired-end $new_fastq_directory/S"$i"_R1_paired.fastq.gz $new_fastq_directory/S"$i"_R2_paired.fastq.gz $output_path/"$plant"_index/"$plant"_index $output_path/calc_expression/$output_suffix"$i"/$output_suffix"$i"
        else
        rsem-calculate-expression -p 20 --bowtie2 $new_fastq_directory/S"$i"_R1.fastq.gz $output_path/"$plant"_index/"$plant"_index $output_path/calc_expression/$output_suffix"$i"/$output_suffix"$i"
        fi

        ### sort bam files
        samtools sort $output_path/calc_expression/$output_suffix"$i"/$output_suffix"$i".transcript.bam -o $output_path/sorted_bam/$output_suffix"$i"_sorted.bam


        ### check if files were created successfully
        if [[ ! -f $output_path/calc_expression/$output_suffix"$i"/$output_suffix"$i".genes.results ]]; then
            echo "Error in 'rsem-calculate-expression' ('$output_suffix"$i".genes.results' file)" >> $output_path/pipeline_$current_date.log
        fi

        if [[ ! -f $output_path/sorted_bam/$output_suffix"$i"_sorted.bam ]]; then
            echo "Error in 'samtools sort' ('$output_suffix"$i"_sorted.bam' file)" >> $output_path/pipeline_$current_date.log
        fi

        echo "calculate expression for "$output_suffix"$i"": done" >> $output_path/pipeline_$current_date.log
done

# copy the '*.genes.results' files to a directory
mkdir $output_path/gene.results.files
cp $output_path/calc_expression/$output_suffixg*/$output_suffix*.genes.results $output_path/gene.results.files

echo "" >> $output_path/pipeline_$current_date.log
date >> $output_path/pipeline_$current_date.log
echo "
▓█████▄  ▒█████   ███▄    █ ▓█████ 
▒██▀ ██▌▒██▒  ██▒ ██ ▀█   █ ▓█   ▀ 
░██   █▌▒██░  ██▒▓██  ▀█ ██▒▒███   
░▓█▄   ▌▒██   ██░▓██▒  ▐▌██▒▒▓█  ▄ 
░▒████▓ ░ ████▓▒░▒██░   ▓██░░▒████▒
 ▒▒▓  ▒ ░ ▒░▒░▒░ ░ ▒░   ▒ ▒ ░░ ▒░ ░
 ░ ▒  ▒   ░ ▒ ▒░ ░ ░░   ░ ▒░ ░ ░  ░
 ░ ░  ░ ░ ░ ░ ▒     ░   ░ ░    ░   
   ░        ░ ░           ░    ░  ░
 ░                                 
"