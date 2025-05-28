#!/bin/bash
usage_yo="
###############################################################################
 YO - 260525
 Bismark WGBS pipeline
 
 ------------------------------------------------------------
  
 Usage:
 ------
 run_bismark_yo.sh [-s <required>] [-g <required>] [options]
 
 Options:
 --------
 -s, --samples       Tab-delimited two-column file: sample-name <TAB> fastq-path
 -g, --genome        FASTA of the reference genome (will be indexed)
 -o, --outdir        Output directory [default: ./bismark_results]
 -n, --ncores        Max CPU cores (multiples of 4 recommended) [default: 16]
     --cx            Produce and keep only '*.CX_report.txt' file
     --sort          Sort & index BAM files (applies only if --cx is off)
     --help

 ------------------------------------------------------------

 Example:
 --------
 Download TAIR10 reference genome:
 ---------------------------------
 cd /PATH/TO
 wget -O TAIR10_chr_all.fas.gz https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
 
 Create a sample table file (example):
 -------------------------------------
 mut_1    PATH/TO/FILE/m1_R1.fastq
 mut_1    PATH/TO/FILE/m1_R2.fastq
 mut_2    PATH/TO/FILE/m2_R1.fastq
 mut_2    PATH/TO/FILE/m2_R2.fastq
 wt_1    PATH/TO/FILE/wt1_R1.fastq
 wt_1    PATH/TO/FILE/wt1_R2.fastq
 wt_2    PATH/TO/FILE/wt2_R1.fastq
 wt_2    PATH/TO/FILE/wt2_R2.fastq
 
 Run:
 ----
 ./run_bismark_yo.sh -s samples_table.txt -g TAIR10_chr_all.fas.gz --cx
###############################################################################
"

### default values
sample_table=
genome_file_full_path=
output_path="./bismark_results"
output_suffix="wgbs_bismark_$(date +%d%m%y)"
n_cores=16
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
        -n | --ncores)
            n_cores=$2
            shift 2
        ;;
        --cx)
            keep_cx=true
            shift
        ;;
        --sort)
            sort_bam=true
            shift
        ;;
        -h | --help)
            echo "$usage_yo"
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

### n-cores for bismark
if [ $n_cores -gt 2 ]; then
    n_cores_2=$((n_cores / 4))
else
    n_cores_2=1
fi

####################

ori_path=$(pwd)
mkdir $output_path
cd $output_path
output_path=$(pwd)

### tmp file for analysis
mkdir $output_path/tmp
cd $output_path/tmp

### Generate log file with a timestamp
log_file="../${output_suffix}.log"
echo "**  $(date +"%d-%m-%y %H:%M")" > "$log_file"
echo "**  samples: ${sample_name[@]}" >> "$log_file"
#echo "paired-end sequence: $paired_end_sequence" >> "$log_file"
if [[ "$keep_cx" == "true" ]]; then
    echo "**  keep just 'CX_report' file" >> "$log_file"
fi
if [[ "$sort_bam" == "true" ]]; then
    echo "**  sort bam file" >> "$log_file"
fi

echo "" >> "$log_file"

### index genom
genome_b_name=$(basename "$genome_file_full_path")
genome_new_path=$output_path/genome_indx/$genome_b_name
echo "indexing genome file: '$genome_b_name'" >> "$log_file"
mkdir -p $output_path/genome_indx
cp $genome_file_full_path $output_path/genome_indx
# Rename genome file if necessary
if [[ "$genome_b_name" == *.fas ]]; then
    mv "$genome_new_path" "${genome_new_path%.fas}.fa"
    echo "* rename genome file: '$(basename ${genome_new_path%.fas}.fa)'" >> "$log_file"
    elif [[ "$genome_new_path" == *.fas.gz ]]; then
    mv "$genome_new_path" "${genome_new_path%.fas.gz}.fa.gz"
    echo "* rename genome file: '$(basename ${genome_new_path%.fas.gz}.fa.gz)'" >> "$log_file"
fi
bismark_genome_preparation $output_path/genome_indx

echo "" >> "$log_file"
echo "-----------------------------------" >> "$log_file"

####################

for ((u = 0; u < ${#sample_name[@]}; u++)); do
    echo "**  $(date +"%d-%m-%y %H:%M")" >> "$log_file"
    echo "" >> "$log_file"

    i="${sample_name[$u]}"
    R1_i="${R1_fastq_path[$u]}"
    R2_i="${R2_fastq_path[$u]}"

    echo "Processing sample: $i" >> "$log_file"
    mkdir -p "$output_path/$i"

    ### mapping to genome
    if [[ "$paired_end_sequence" == "false" ]]; then
        Rs_type="se"
        echo "mapping to genome for single-end sequence:" >> "$log_file"
        echo "read1 file: $R1_i" >> "$log_file"
        bismark --bowtie2 --parallel "$n_cores_2" $output_path/genome_indx "$R1_i" -o $output_path/"$i" --prefix "$i" # --basename
    else
        Rs_type="pe"
        echo "mapping to genome for peired-end sequence:" >> "$log_file"
        echo "* read1 file: '$(basename "$R1_i")'" >> "$log_file"
        echo "* read2 file: '$(basename "$R2_i")'" >> "$log_file"
        bismark --bowtie2 --parallel "$n_cores_2" $output_path/genome_indx -1 "$R1_i" -2 "$R2_i" -o $output_path/"$i" --prefix "$i" # --basename
    fi
    mv $output_path/"$i"/"$i"*.bam $output_path/"$i"/"$i"_bismark_"$Rs_type".bam # change name
    mv $output_path/"$i"/"$i"*_report.txt $output_path/"$i"/"$i"_bismark_"$Rs_type"_report.txt # change name

    ### methylation calling
    echo "" >> "$log_file"
    echo "methylation calling..." >> "$log_file"
    mkdir -p $output_path/"$i"/methylation_extractor

    if [[ "$keep_cx" == "true" ]]; then
        # run 'methylation_extractor' and keep 'CX_report' file only
        bismark_methylation_extractor --cytosine_report --CX --parallel "$n_cores_2" --genome_folder $output_path/genome_indx -o $output_path/"$i"/methylation_extractor $output_path/"$i"/"$i"_bismark_"$Rs_type".bam

        mv  $output_path/"$i"/methylation_extractor/*.CX_report.txt $output_path

        # check if CX_report file exists
        if [[ -f "$output_path/"$i"/methylation_extractor/*.CX_report.txt" ]]; then
            echo "CX_report file: '"$i"_bismark_"$Rs_type".CX_report.txt'" >> "$log_file"
        else
            echo "Error: CX_report file: '"$i"_bismark_"$Rs_type".CX_report.txt' does not exist." >> "$log_file"
        fi

        rm -r $output_path/"$i"

    else
        # run 'methylation_extractor' and keep all files (without 'CX_report')
        bismark_methylation_extractor --parallel "$n_cores_2" --genome_folder $output_path/genome_indx -o $output_path/"$i"/methylation_extractor $output_path/"$i"/"$i"_bismark_"$Rs_type".bam
        echo "output files: '$output_path/"$i"/methylation_extractor/'" >> "$log_file"

        # sort bam files (can use in IGV softwar)
        if [[ "$sort_bam" == "true" ]]; then
            samtools sort $output_path/"$i"/"$i"_*.bam -o $output_path/"$i"/"$i"_*_sorted.bam
            samtools index $output_path/"$i"/"$i"_*_sorted.bam
        fi
    fi

    echo "" >> "$log_file"
    echo "-----------------------------------" >> "$log_file"
done

if [[ "$keep_cx" == "true" ]]; then
    rm -r $output_path/genome_indx
fi

echo "**  $(date +"%d-%m-%y %H:%M")" >> "$log_file"
cd $ori_path
rm -r $output_path/tmp
