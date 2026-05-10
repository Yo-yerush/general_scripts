#!/bin/bash

usage_yo="
###############################################################################
YO - 091025
RSEM RNA-seq pipeline (Bowtie2)

------------------------------------------------------------

Usage:
------
run_rsem_yo.sh [-s <required>] [-g <required>] [-a <required>] [options]

Options:
--------
-s, --samples       Tab-delimited two-column file: sample-name <TAB> fastq-path
-g, --genome        FASTA of the reference genome (will be indexed)
-a, --gtf           GTF/GFF3 annotation file (will be copied next to the index)
-o, --outdir        Output directory [default: ./rsem_results]
-n, --ncores        Number of maximum cores usage for alignment/quant [default: 16]
-c, --chunks        How many samples to run in parallel [default: auto]
--sort              Sort & index transcript-aligned BAM files
--genes_results     Keep just '*.genes.results' file (for downstream analysis)
--help

Notes:
------
* Sample table supports single- or paired-end reads. Use typical names like:
_R1_, _R2_, _1.fq(.gz), _2.fastq(.gz). If single-end, just list one row per
sample with the read file path.
* Keeps the overall UX and logging style of your Bismark pipeline, but shorter.

Example:
--------
Create a sample table file (example):
-------------------------------------
mt_1    /PATH/TO/mt1_R1.fastq.gz
mt_1    /PATH/TO/mt1_R2.fastq.gz
wt_1    /PATH/TO/wt1_R1.fastq.gz
wt_1    /PATH/TO/wt1_R2.fastq.gz

Run:
----
$ ./run_rsem_yo.sh -s samples_table.txt -g TAIR10_chr_all.fa -a TAIR10.gtf --genes_results
###############################################################################
"

####################
### default values
sample_table=
genome_file_full_path=
ann_file_full_path=
output_path="./rsem_results"
output_suffix="rnaseq_rsem_$(date +%d%m%y)"
n_cores=16
chunks=auto
sort_bam=false
genes_res=false

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
        -a | --gtf)
            ann_file_full_path=$2
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
        -c | --chunks)
            chunks=$2
            shift 2
        ;;
        --sort)
            sort_bam=true
            shift
        ;;
        --genes_results)
            genes_res=true
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
# basic checks
if [[ -z "$sample_table" || -z "$genome_file_full_path" || -z "$ann_file_full_path" ]]; then
    echo "Error: -s/--samples, -g/--genome, and -a/--gtf are required."
    echo ""
    echo "$usage_yo"
    exit 1
fi

if [[ ! -f "$sample_table" ]]; then
    echo "Error: Sample table file '$sample_table' does not exist."
    exit 1
fi
if [[ ! -f "$genome_file_full_path" ]]; then
    echo "Error: Genome file '$genome_file_full_path' does not exist."
    exit 1
fi
if [[ ! -f "$ann_file_full_path" ]]; then
    echo "Error: Annotation (GTF/GFF3) file '$ann_file_full_path' does not exist."
    exit 1
fi

if [[ "$sort_bam" == "true" && "$genes_res" == "true" ]]; then
    echo "Error: --sort and --genes_results options cannot be used together."
    echo "Sorting BAM files requires keeping the full output directory structure."
    exit 1
fi

# ensure sample table has unix line endings
dos2unix "$sample_table" 2>/dev/null

####################
### read sample names and fastq file paths as an array
mapfile -t sample_name < <(awk '!seen[$1]++ {print $1}' "$sample_table")
mapfile -t R1_fastq_path < <(awk '$2 ~ /(_R1_|_R1\.fq|_R1\.fastq|_1\.fq|_1\.fastq)/ {print $2}' "$sample_table")
mapfile -t R2_fastq_path < <(awk '$2 ~ /(_R2_|_R2\.fq|_R2\.fastq|_2\.fq|_2\.fastq)/ {print $2}' "$sample_table")

paired_end_sequence=true
((${#R2_fastq_path[@]} == 0)) && paired_end_sequence=false
if [[ "$paired_end_sequence" == "false" ]]; then
    # For single-end, use all rows in the samples table as read file paths
    mapfile -t R1_fastq_path < <(awk '{print $2}' "$sample_table")
fi

####################
### chunking and cores per job
sample_count=${#sample_name[@]}

if [[ "$chunks" == "auto" ]]; then
    if (( sample_count >= 12 )); then ### and module!!
        chunks=4
    elif (( sample_count >= 6 )); then
        chunks=3
    else
        chunks=2
    fi
fi

# cap chunks to sample count and enforce minimum 1
if (( chunks > sample_count )); then
    chunks=$sample_count
fi
if (( chunks < 1 )); then
    chunks=1
fi

cores_per_job=$(( n_cores / chunks ))
if (( cores_per_job < 1 )); then
    cores_per_job=1
fi

####################
### prepare output
ori_path=$(pwd)
mkdir -p "$output_path"
cd "$output_path"
output_path=$(pwd)

### tmp folder + log
mkdir -p "$output_path/tmp"
cd "$output_path/tmp"

echo ""
echo "**  samples: ${sample_name[@]}"
if [[ "$sort_bam" == "true" ]]; then
    echo "**  will sort & index transcript BAMs"
fi
if [[ "$genes_res" == "true" ]]; then
    echo "**  will keep just '*.genes.results' files"
fi
echo ""

mkdir -p "${output_path}/logs"

####################
### index (prepare) the reference for RSEM
mkdir -p "$output_path/genome_indx"
genome_b_name=$(basename "$genome_file_full_path")
ann_b_name=$(basename "$ann_file_full_path")

# copy inputs
cp "$genome_file_full_path" "$output_path/genome_indx"/
cp "$ann_file_full_path" "$output_path/genome_indx"/

genome_new_path="$output_path/genome_indx/$genome_b_name"
ann_new_path="$output_path/genome_indx/$ann_b_name"

# decompress if gz
if [[ "$genome_new_path" == *.gz ]]; then
    gunzip -c "$genome_new_path" > "${genome_new_path%.gz}"
    genome_new_path="${genome_new_path%.gz}"
fi
if [[ "$ann_new_path" == *.gz ]]; then
    gunzip -c "$ann_new_path" > "${ann_new_path%.gz}"
    ann_new_path="${ann_new_path%.gz}"
fi

# rename .fas -> .fa (cosmetic) and handle gz
if [[ "$genome_new_path" == *.fas ]]; then
    mv "$genome_new_path" "${genome_new_path%.fas}.fa"
    genome_new_path="${genome_new_path%.fas}.fa"
    echo "** rename genome file: '$(basename "$genome_new_path")'"
fi

index_prefix="$output_path/genome_indx/rsem_index"
cd "$output_path/genome_indx"
echo "preparing RSEM reference..."
echo "**  $(date +"%d-%m-%y %H:%M")" > "${output_path}/logs/index.rsem.log"
rsem-prepare-reference --bowtie2 --gtf "$(basename "$ann_new_path")" "$(basename "$genome_new_path")" "rsem_index" >> "${output_path}/logs/index.rsem.log" 2>&1

cd "$output_path/tmp"
echo ""
echo "-----------------------------------"

####################
### main loop
for ((u = 0; u < ${#sample_name[@]}; u++)); do
    echo "**  $(date +"%d-%m-%y %H:%M")"
    echo ""

    i="${sample_name[$u]}"
    R1_i="${R1_fastq_path[$u]}"
    R2_i="${R2_fastq_path[$u]}"

    echo "Processing sample: $i"
    mkdir -p "$output_path/$i"

    if [[ "$paired_end_sequence" == "false" ]]; then
        echo "* read1 file: '$(basename "$R1_i")'"
    else
        echo "* read1 file: '$(basename "$R1_i")'"
        echo "* read2 file: '$(basename "$R2_i")'"
    fi

    echo "**  $(date +"%d-%m-%y %H:%M")" > "${output_path}/logs/sample_${i}.rsem.log"

    (
        echo ""
        
        if [[ "$paired_end_sequence" == "false" ]]; then
            rsem-calculate-expression -p "$cores_per_job" --bowtie2 "$R1_i" "$index_prefix" "$output_path/$i/$i"
        else
            rsem-calculate-expression -p "$cores_per_job" --bowtie2 --paired-end "$R1_i" "$R2_i" "$index_prefix" "$output_path/$i/$i"
        fi

        if [[ "$sort_bam" == "true" ]]; then
            if [[ -f "$output_path/$i/$i.transcript.bam" ]]; then
                samtools sort "$output_path/$i/$i.transcript.bam" -o "$output_path/$i/${i}_sorted.bam"
                samtools index "$output_path/$i/${i}_sorted.bam"
            fi
        fi

        if [[ "$genes_res" == "true" ]]; then
            mv "$output_path/$i/$i".genes.results "$output_path"
            rm -r -- "$output_path/$i"
        fi
        
        echo ""
        echo "**  $(date +"%d-%m-%y %H:%M")"
    ) >> "${output_path}/logs/sample_${i}.rsem.log" 2>&1 &


if (( (u + 1) % chunks == 0 )); then
    wait
fi

    echo ""
    echo "-----------------------------------"
done

wait

echo "**  $(date +"%d-%m-%y %H:%M")"
cd "$ori_path"
rm -r "$output_path/tmp"