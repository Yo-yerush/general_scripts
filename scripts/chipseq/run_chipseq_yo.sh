#!/usr/bin/env bash

set -euo pipefail

usage_yo="
###############################################################################
YO - 300726
Arabidopsis ChIP-seq pipeline (Bowtie2 + MACS3 + DiffBind)

------------------------------------------------------------

Usage:
------
run_chipseq_yo.sh -s <samples.tsv> -g <genome.fa[.gz]> \\
  --contrast <TEST,REFERENCE> [options]

Required:
---------
-s, --samples           Tab-delimited sample table with this exact header/order:
                        SampleID  Condition  Replicate  Factor  Type  ControlID  Read1  Read2
-g, --genome            Reference genome FASTA; it will be copied, decompressed,
                        and indexed with Bowtie2.
--contrast              Differential comparison as TEST,REFERENCE.
                        Example: --contrast mto3,WT
                        Positive DiffBind Fold values mean higher enrichment in TEST.

Sample-table rules:
-------------------
* One row per sequencing library.
* Type must be ChIP or Input.
* ControlID in each ChIP row must match an Input SampleID.
* Use '.' in Read2 for single-end data.
* Analyze one antibody/histone mark per run (one Factor among ChIP rows).
* Differential analysis requires at least two biological ChIP replicates per condition.

Options:
--------
-o, --outdir            Output directory [default: ./chipseq_results]
-n, --ncores            Maximum total cores [default: 16]
-c, --chunks            Samples processed in parallel [default: auto]
--trim                   Adapter/quality trimming with fastp [default: off]
--fastqc                 Run FastQC on reads used for mapping [default: off]
--broad                  Call broad peaks (e.g. H3K9me2/H3K27me3)
--qvalue                 MACS3 narrow-peak q-value [default: 0.01]
--broad-cutoff           MACS3 broad-region q-value [default: 0.10]
--gsize                  Effective genome size for MACS3/RPGC; integer or auto
                        [default: auto = count non-N bases in FASTA]
--mapq                   Minimum mapping quality [default: 30]
--max-insert             Maximum valid paired-end fragment length [default: 2000]
--exclude-contigs        Comma-separated contigs excluded after mapping
                        [default: ChrC,ChrM; use 'none' to keep all]
--keep-duplicates        Retain PCR/optical duplicates [default: remove]
--normalization          bigWig normalization: CPM or RPGC [default: CPM]
--bin-size               bigWig bin size [default: 10]
--min-overlap            Minimum peaksets supporting a DiffBind consensus site
                        [default: 2]
--fdr                    FDR threshold for significant differential peaks
                        [default: 0.05]
--lfc                    Minimum absolute log2 fold-change for the significant table
                        [default: 0]
--no-diff                Stop after mapping, tracks, and peak calling
--keep-intermediate      Keep pre-filter and duplicate-marked BAM files
-h, --help

Example sample table:
---------------------
SampleID\tCondition\tReplicate\tFactor\tType\tControlID\tRead1\tRead2
WT_H3K9me2_1\tWT\t1\tH3K9me2\tChIP\tWT_Input_1\t/path/WT_K9_1_R1.fastq.gz\t/path/WT_K9_1_R2.fastq.gz
WT_Input_1\tWT\t1\tInput\tInput\t.\t/path/WT_Input_1_R1.fastq.gz\t/path/WT_Input_1_R2.fastq.gz
mto3_H3K9me2_1\tmto3\t1\tH3K9me2\tChIP\tmto3_Input_1\t/path/mto3_K9_1_R1.fastq.gz\t/path/mto3_K9_1_R2.fastq.gz
mto3_Input_1\tmto3\t1\tInput\tInput\t.\t/path/mto3_Input_1_R1.fastq.gz\t/path/mto3_Input_1_R2.fastq.gz

Example run for a broad histone mark:
-------------------------------------
./run_chipseq_yo.sh -s chipseq_samples.tsv -g TAIR10_chr_all.fa.gz \\
  --contrast mto3,WT --broad --trim -n 16 -o H3K9me2_results

Main software:
--------------
bowtie2, samtools, picard, macs3, deepTools (bamCoverage/bamCompare), R/DiffBind
Optional: fastp (--trim), fastqc (--fastqc), multiqc (run automatically if installed)
###############################################################################
"

####################
### default values
sample_table=""
genome_file_full_path=""
output_path="./chipseq_results"
contrast_arg=""
n_cores=16
chunks="auto"
run_trim=false
run_fastqc=false
broad_peaks=false
macs_qvalue="0.01"
broad_cutoff="0.10"
effective_genome_size="auto"
mapq=30
max_insert=2000
exclude_contigs="ChrC,ChrM"
keep_duplicates=false
normalization="CPM"
bin_size=10
min_overlap=2
fdr_threshold="0.05"
lfc_threshold="0"
run_diff=true
keep_intermediate=false

while [[ $# -gt 0 ]]; do
    case "$1" in
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
        -c | --chunks)
            chunks=$2
            shift 2
        ;;
        --contrast)
            contrast_arg=$2
            shift 2
        ;;
        --trim)
            run_trim=true
            shift
        ;;
        --fastqc)
            run_fastqc=true
            shift
        ;;
        --broad)
            broad_peaks=true
            shift
        ;;
        --qvalue)
            macs_qvalue=$2
            shift 2
        ;;
        --broad-cutoff)
            broad_cutoff=$2
            shift 2
        ;;
        --gsize)
            effective_genome_size=$2
            shift 2
        ;;
        --mapq)
            mapq=$2
            shift 2
        ;;
        --max-insert)
            max_insert=$2
            shift 2
        ;;
        --exclude-contigs)
            exclude_contigs=$2
            shift 2
        ;;
        --keep-duplicates)
            keep_duplicates=true
            shift
        ;;
        --normalization)
            normalization=$2
            shift 2
        ;;
        --bin-size)
            bin_size=$2
            shift 2
        ;;
        --min-overlap)
            min_overlap=$2
            shift 2
        ;;
        --fdr)
            fdr_threshold=$2
            shift 2
        ;;
        --lfc)
            lfc_threshold=$2
            shift 2
        ;;
        --no-diff)
            run_diff=false
            shift
        ;;
        --keep-intermediate)
            keep_intermediate=true
            shift
        ;;
        -h | --help)
            echo "$usage_yo"
            exit 0
        ;;
        *)
            echo "Error: Unknown option: $1" >&2
            echo "$usage_yo" >&2
            exit 1
        ;;
    esac
done

####################
### helper functions
fail() {
    echo "Error: $*" >&2
    exit 1
}

log() {
    echo "**  $(date +"%d-%m-%y %H:%M") | $*"
}

command_required() {
    command -v "$1" >/dev/null 2>&1 || fail "Required command '$1' was not found in PATH."
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

absolute_file() {
    local path=$1
    if command -v realpath >/dev/null 2>&1; then
        realpath "$path"
    else
        readlink -f "$path"
    fi
}

wait_for_batch() {
    local -n batch_pids_ref=$1
    local -n batch_labels_ref=$2
    local failed=0
    local idx

    for idx in "${!batch_pids_ref[@]}"; do
        if ! wait "${batch_pids_ref[$idx]}"; then
            echo "Error: Job failed for '${batch_labels_ref[$idx]}'. See its log file." >&2
            failed=1
        fi
    done
    (( failed == 0 )) || fail "One or more parallel jobs failed."
}

####################
### basic argument checks
[[ -n "$sample_table" ]] || fail "-s/--samples is required."
[[ -n "$genome_file_full_path" ]] || fail "-g/--genome is required."
[[ -f "$sample_table" ]] || fail "Sample table '$sample_table' does not exist."
[[ -f "$genome_file_full_path" ]] || fail "Genome FASTA '$genome_file_full_path' does not exist."

if [[ "$run_diff" == "true" && -z "$contrast_arg" ]]; then
    fail "--contrast TEST,REFERENCE is required unless --no-diff is used."
fi

is_positive_integer "$n_cores" || fail "--ncores must be a positive integer."
[[ "$chunks" == "auto" ]] || is_positive_integer "$chunks" || fail "--chunks must be 'auto' or a positive integer."
is_positive_integer "$mapq" || [[ "$mapq" == "0" ]] || fail "--mapq must be a non-negative integer."
is_positive_integer "$max_insert" || fail "--max-insert must be a positive integer."
is_positive_integer "$bin_size" || fail "--bin-size must be a positive integer."
is_positive_integer "$min_overlap" || fail "--min-overlap must be a positive integer."

normalization=$(echo "$normalization" | tr '[:lower:]' '[:upper:]')
[[ "$normalization" == "CPM" || "$normalization" == "RPGC" ]] || fail "--normalization must be CPM or RPGC."

if [[ "$effective_genome_size" != "auto" ]]; then
    is_positive_integer "$effective_genome_size" || fail "--gsize must be 'auto' or a positive integer."
fi

if [[ "$run_diff" == "true" ]]; then
    [[ "$contrast_arg" == *,* ]] || fail "--contrast must be formatted as TEST,REFERENCE."
    test_condition=${contrast_arg%%,*}
    reference_condition=${contrast_arg#*,}
    [[ -n "$test_condition" && -n "$reference_condition" ]] || fail "Both TEST and REFERENCE must be supplied to --contrast."
    [[ "$test_condition" != "$reference_condition" ]] || fail "TEST and REFERENCE conditions must differ."
else
    test_condition=""
    reference_condition=""
fi

####################
### required software
for cmd in bowtie2 bowtie2-build samtools picard macs3 bamCoverage bamCompare awk sed grep sort uniq gzip Rscript; do
    command_required "$cmd"
done
[[ "$run_trim" == "false" ]] || command_required fastp
[[ "$run_fastqc" == "false" ]] || command_required fastqc

if [[ "$run_diff" == "true" ]]; then
    Rscript -e 'if (!requireNamespace("DiffBind", quietly=TRUE)) quit(status=2)' \
        >/dev/null 2>&1 || fail "R package 'DiffBind' is not installed. Install it with BiocManager::install('DiffBind')."
fi

####################
### normalize input paths and line endings
sample_table=$(absolute_file "$sample_table")
genome_file_full_path=$(absolute_file "$genome_file_full_path")
original_sample_table="$sample_table"

# Work from a clean temporary copy; do not modify the user's table.
clean_sample_table=$(mktemp)
sed 's/\r$//' "$sample_table" > "$clean_sample_table"
sample_table="$clean_sample_table"
trap 'rm -f "$clean_sample_table"' EXIT

expected_header=$'SampleID\tCondition\tReplicate\tFactor\tType\tControlID\tRead1\tRead2'
actual_header=$(head -n 1 "$sample_table")
[[ "$actual_header" == "$expected_header" ]] || fail "Sample-table header/order is invalid. Expected: $expected_header"

####################
### read sample table
mapfile -t sample_id < <(awk -F '\t' 'NR>1 && NF>0 && $1 !~ /^#/ {print $1}' "$sample_table")
mapfile -t condition < <(awk -F '\t' 'NR>1 && NF>0 && $1 !~ /^#/ {print $2}' "$sample_table")
mapfile -t replicate < <(awk -F '\t' 'NR>1 && NF>0 && $1 !~ /^#/ {print $3}' "$sample_table")
mapfile -t factor < <(awk -F '\t' 'NR>1 && NF>0 && $1 !~ /^#/ {print $4}' "$sample_table")
mapfile -t library_type < <(awk -F '\t' 'NR>1 && NF>0 && $1 !~ /^#/ {print $5}' "$sample_table")
mapfile -t control_id < <(awk -F '\t' 'NR>1 && NF>0 && $1 !~ /^#/ {print $6}' "$sample_table")
mapfile -t read1_path < <(awk -F '\t' 'NR>1 && NF>0 && $1 !~ /^#/ {print $7}' "$sample_table")
mapfile -t read2_path < <(awk -F '\t' 'NR>1 && NF>0 && $1 !~ /^#/ {print $8}' "$sample_table")

sample_count=${#sample_id[@]}
(( sample_count > 0 )) || fail "No sample rows were found."

# Require eight populated fields; Read2 may be '.'.
invalid_nf=$(awk -F '\t' 'NR>1 && NF>0 && $1 !~ /^#/ && NF!=8 {print NR}' "$sample_table" | paste -sd, -)
[[ -z "$invalid_nf" ]] || fail "Sample table must contain exactly 8 tab-delimited columns. Invalid line(s): $invalid_nf"

duplicate_ids=$(printf '%s\n' "${sample_id[@]}" | sort | uniq -d | paste -sd, -)
[[ -z "$duplicate_ids" ]] || fail "Duplicate SampleID values: $duplicate_ids"

# Associative arrays allow ChIP rows to find their matching Input library.
declare -A type_by_id
declare -A condition_by_id
declare -A replicate_by_id
declare -A factor_by_id
declare -A control_by_id
declare -A r1_by_id
declare -A r2_by_id
declare -A paired_by_id

for ((u=0; u<sample_count; u++)); do
    sid=${sample_id[$u]}
    cond=${condition[$u]}
    rep=${replicate[$u]}
    fac=${factor[$u]}
    typ=$(echo "${library_type[$u]}" | tr '[:lower:]' '[:upper:]')
    ctrl=${control_id[$u]}
    r1=${read1_path[$u]}
    r2=${read2_path[$u]}

    [[ -n "$sid" && -n "$cond" && -n "$rep" && -n "$fac" && -n "$typ" && -n "$ctrl" && -n "$r1" ]] \
        || fail "Missing required value in sample row for SampleID '$sid'."
    [[ "$typ" == "CHIP" || "$typ" == "INPUT" ]] || fail "Type for '$sid' must be ChIP or Input."
    [[ -f "$r1" ]] || fail "Read1 for '$sid' does not exist: $r1"
    r1=$(absolute_file "$r1")

    paired=false
    if [[ -n "$r2" && "$r2" != "." && "$r2" != "NA" ]]; then
        [[ -f "$r2" ]] || fail "Read2 for '$sid' does not exist: $r2"
        r2=$(absolute_file "$r2")
        paired=true
    else
        r2="."
    fi

    type_by_id[$sid]=$typ
    condition_by_id[$sid]=$cond
    replicate_by_id[$sid]=$rep
    factor_by_id[$sid]=$fac
    control_by_id[$sid]=$ctrl
    r1_by_id[$sid]=$r1
    r2_by_id[$sid]=$r2
    paired_by_id[$sid]=$paired
done

# Validate control references and collect ChIP libraries.
chip_ids=()
for sid in "${sample_id[@]}"; do
    if [[ "${type_by_id[$sid]}" == "CHIP" ]]; then
        chip_ids+=("$sid")
        ctrl=${control_by_id[$sid]}
        [[ -n "${type_by_id[$ctrl]+x}" ]] || fail "ControlID '$ctrl' for ChIP sample '$sid' is absent from the sample table."
        [[ "${type_by_id[$ctrl]}" == "INPUT" ]] || fail "ControlID '$ctrl' for '$sid' must refer to a Type=Input row."
        [[ "${paired_by_id[$sid]}" == "${paired_by_id[$ctrl]}" ]] || fail "ChIP '$sid' and Input '$ctrl' must both be single-end or both paired-end."
    fi
done
(( ${#chip_ids[@]} > 0 )) || fail "At least one Type=ChIP row is required."

if [[ "$run_diff" == "true" ]]; then
    # Compare one factor only, preventing accidental pooling of different antibodies.
    mapfile -t compared_factors < <(
        for sid in "${chip_ids[@]}"; do
            cond=${condition_by_id[$sid]}
            if [[ "$cond" == "$test_condition" || "$cond" == "$reference_condition" ]]; then
                echo "${factor_by_id[$sid]}"
            fi
        done | sort -u
    )
    (( ${#compared_factors[@]} == 1 )) || fail "The selected contrast must contain exactly one ChIP Factor. Found: ${compared_factors[*]:-none}"

    test_n=0
    reference_n=0
    for sid in "${chip_ids[@]}"; do
        [[ "${condition_by_id[$sid]}" == "$test_condition" ]] && ((test_n+=1))
        [[ "${condition_by_id[$sid]}" == "$reference_condition" ]] && ((reference_n+=1))
    done
    (( test_n >= 2 )) || fail "Condition '$test_condition' has $test_n ChIP replicate(s); DiffBind requires at least 2."
    (( reference_n >= 2 )) || fail "Condition '$reference_condition' has $reference_n ChIP replicate(s); DiffBind requires at least 2."
fi

####################
### chunking and cores per sample
if [[ "$chunks" == "auto" ]]; then
    if (( sample_count >= 12 )); then
        chunks=4
    elif (( sample_count >= 6 )); then
        chunks=3
    else
        chunks=2
    fi
fi
(( chunks > sample_count )) && chunks=$sample_count
(( chunks < 1 )) && chunks=1
cores_per_job=$(( n_cores / chunks ))
(( cores_per_job < 1 )) && cores_per_job=1

####################
### output folders and run log
mkdir -p "$output_path"
output_path=$(cd "$output_path" && pwd)
mkdir -p "$output_path"/{genome_index,alignments,tracks,peaks,qc,logs,tmp,diffbind}
cp "$sample_table" "$output_path/input_samples.tsv"
run_log="$output_path/logs/run_chipseq_yo.log"
exec > >(tee -a "$run_log") 2>&1

log "Starting ChIP-seq pipeline"
echo "**  samples: ${sample_id[*]}"
echo "**  ChIP samples: ${chip_ids[*]}"
echo "**  cores: $n_cores total; $chunks parallel job(s); $cores_per_job core(s)/job"
echo "**  peak mode: $([[ "$broad_peaks" == "true" ]] && echo broad || echo narrow)"
echo "**  duplicate handling: $([[ "$keep_duplicates" == "true" ]] && echo retain || echo remove)"
echo "**  excluded contigs: $exclude_contigs"
echo "**  bigWig normalization: $normalization"
if [[ "$run_diff" == "true" ]]; then
    echo "**  differential contrast: $test_condition vs $reference_condition"
fi
echo ""

####################
### software versions
{
    echo "Run date: $(date -Is)"
    bowtie2 --version | head -n 1
    samtools --version | head -n 1
    picard MarkDuplicates --version 2>&1 | head -n 1 || true
    macs3 --version
    bamCoverage --version
    Rscript --version 2>&1
    if [[ "$run_diff" == "true" ]]; then
        Rscript -e 'cat("DiffBind ", as.character(packageVersion("DiffBind")), "\n", sep="")'
    fi
    [[ "$run_trim" == "false" ]] || fastp --version 2>&1 | head -n 1
    [[ "$run_fastqc" == "false" ]] || fastqc --version 2>&1 | head -n 1
} > "$output_path/logs/software_versions.txt"

####################
### prepare and index genome
genome_fa="$output_path/genome_index/reference.fa"
if [[ "$genome_file_full_path" == *.gz ]]; then
    log "Decompressing reference genome"
    gzip -dc "$genome_file_full_path" > "$genome_fa"
else
    cp "$genome_file_full_path" "$genome_fa"
fi
[[ -s "$genome_fa" ]] || fail "Prepared genome FASTA is empty."

if [[ "$effective_genome_size" == "auto" ]]; then
    effective_genome_size=$(awk '
        /^>/ {next}
        {
            gsub(/[Nn[:space:]]/, "", $0)
            total += length($0)
        }
        END {print total+0}
    ' "$genome_fa")
    (( effective_genome_size > 0 )) || fail "Could not calculate effective genome size."
fi
echo "$effective_genome_size" > "$output_path/genome_index/effective_genome_size.txt"
echo "**  effective genome size: $effective_genome_size"

index_prefix="$output_path/genome_index/bowtie2_index"
log "Building Bowtie2 index"
bowtie2-build --threads "$n_cores" "$genome_fa" "$index_prefix" \
    > "$output_path/logs/bowtie2_build.log" 2>&1
samtools faidx "$genome_fa"

####################
### process one sequencing library
process_library() {
    local sid=$1
    local r1=${r1_by_id[$sid]}
    local r2=${r2_by_id[$sid]}
    local paired=${paired_by_id[$sid]}
    local sample_dir="$output_path/alignments/$sid"
    local sample_qc="$output_path/qc/$sid"
    local sample_tmp="$output_path/tmp/$sid"
    local map_log="$output_path/logs/${sid}.mapping.log"
    local r1_used=$r1
    local r2_used=$r2

    mkdir -p "$sample_dir" "$sample_qc" "$sample_tmp"
    echo "**  $(date +"%d-%m-%y %H:%M") | Processing library: $sid"

    if [[ "$run_trim" == "true" ]]; then
        if [[ "$paired" == "true" ]]; then
            r1_used="$sample_tmp/${sid}_R1.trimmed.fastq.gz"
            r2_used="$sample_tmp/${sid}_R2.trimmed.fastq.gz"
            fastp \
                -i "$r1" -I "$r2" \
                -o "$r1_used" -O "$r2_used" \
                --detect_adapter_for_pe \
                --thread "$cores_per_job" \
                --html "$sample_qc/${sid}.fastp.html" \
                --json "$sample_qc/${sid}.fastp.json" \
                > "$output_path/logs/${sid}.fastp.log" 2>&1
        else
            r1_used="$sample_tmp/${sid}.trimmed.fastq.gz"
            fastp \
                -i "$r1" -o "$r1_used" \
                --thread "$cores_per_job" \
                --html "$sample_qc/${sid}.fastp.html" \
                --json "$sample_qc/${sid}.fastp.json" \
                > "$output_path/logs/${sid}.fastp.log" 2>&1
        fi
    fi

    if [[ "$run_fastqc" == "true" ]]; then
        if [[ "$paired" == "true" ]]; then
            fastqc -t "$cores_per_job" -o "$sample_qc" "$r1_used" "$r2_used" \
                > "$output_path/logs/${sid}.fastqc.log" 2>&1
        else
            fastqc -t "$cores_per_job" -o "$sample_qc" "$r1_used" \
                > "$output_path/logs/${sid}.fastqc.log" 2>&1
        fi
    fi

    local sorted_bam="$sample_dir/${sid}.sorted.bam"
    if [[ "$paired" == "true" ]]; then
        bowtie2 --very-sensitive --no-mixed --no-discordant -X "$max_insert" \
            -x "$index_prefix" -1 "$r1_used" -2 "$r2_used" -p "$cores_per_job" \
            2> "$map_log" \
            | samtools view -@ "$cores_per_job" -b - \
            | samtools sort -@ "$cores_per_job" -o "$sorted_bam" -
    else
        bowtie2 --very-sensitive \
            -x "$index_prefix" -U "$r1_used" -p "$cores_per_job" \
            2> "$map_log" \
            | samtools view -@ "$cores_per_job" -b - \
            | samtools sort -@ "$cores_per_job" -o "$sorted_bam" -
    fi
    samtools index -@ "$cores_per_job" "$sorted_bam"
    samtools flagstat -@ "$cores_per_job" "$sorted_bam" > "$sample_qc/${sid}.before_filter.flagstat.txt"

    local marked_bam="$sample_dir/${sid}.marked_duplicates.bam"
    mkdir -p "$sample_tmp/picard"
    picard MarkDuplicates \
        I="$sorted_bam" \
        O="$marked_bam" \
        M="$sample_qc/${sid}.duplicate_metrics.txt" \
        CREATE_INDEX=true \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR="$sample_tmp/picard" \
        > "$output_path/logs/${sid}.picard.log" 2>&1

    local exclude_flags=1796  # unmapped, secondary, QC-fail, duplicates
    [[ "$keep_duplicates" == "false" ]] || exclude_flags=772

    local final_bam="$sample_dir/${sid}.filtered.bam"
    local proper_pair_args=()
    [[ "$paired" == "false" ]] || proper_pair_args=(-f 2)

    local region_args=()
    if [[ -n "$exclude_contigs" && "$exclude_contigs" != "none" && "$exclude_contigs" != "NONE" ]]; then
        local exclude_regex
        exclude_regex="^($(echo "$exclude_contigs" | sed 's/,/|/g'))$"
        mapfile -t keep_contigs < <(
            samtools idxstats "$marked_bam" \
                | awk -v re="$exclude_regex" '$1!="*" && $1 !~ re {print $1}'
        )
        (( ${#keep_contigs[@]} > 0 )) || fail "No contigs remained after excluding '$exclude_contigs'."
        region_args=("${keep_contigs[@]}")
    fi

    samtools view -@ "$cores_per_job" -b -q "$mapq" \
        "${proper_pair_args[@]}" -F "$exclude_flags" \
        -o "$final_bam" "$marked_bam" "${region_args[@]}"
    samtools index -@ "$cores_per_job" "$final_bam"
    samtools flagstat -@ "$cores_per_job" "$final_bam" > "$sample_qc/${sid}.filtered.flagstat.txt"
    samtools idxstats "$final_bam" > "$sample_qc/${sid}.filtered.idxstats.txt"

    local coverage_args=(
        -b "$final_bam"
        -o "$output_path/tracks/${sid}.${normalization}.bw"
        --normalizeUsing "$normalization"
        --binSize "$bin_size"
        --numberOfProcessors "$cores_per_job"
    )
    if [[ "$normalization" == "RPGC" ]]; then
        coverage_args+=(--effectiveGenomeSize "$effective_genome_size")
    fi
    if [[ "$paired" == "false" ]]; then
        coverage_args+=(--extendReads 200)
    fi
    bamCoverage "${coverage_args[@]}" > "$output_path/logs/${sid}.bamCoverage.log" 2>&1

    if [[ "$keep_intermediate" == "false" ]]; then
        rm -f "$sorted_bam" "$sorted_bam.bai" "$marked_bam" "$marked_bam.bai"
        rm -rf "$sample_tmp"
    fi

    echo "**  $(date +"%d-%m-%y %H:%M") | Finished library: $sid"
}

####################
### mapping loop
log "Mapping and filtering libraries"
batch_pids=()
batch_labels=()
for ((u=0; u<sample_count; u++)); do
    sid=${sample_id[$u]}
    process_library "$sid" > "$output_path/logs/${sid}.pipeline.log" 2>&1 &
    batch_pids+=("$!")
    batch_labels+=("$sid")

    if (( ${#batch_pids[@]} == chunks )); then
        wait_for_batch batch_pids batch_labels
        batch_pids=()
        batch_labels=()
    fi
done
if (( ${#batch_pids[@]} > 0 )); then
    wait_for_batch batch_pids batch_labels
fi

if command -v multiqc >/dev/null 2>&1; then
    log "Running MultiQC"
    multiqc -f -o "$output_path/qc/multiqc" "$output_path/qc" "$output_path/logs" \
        > "$output_path/logs/multiqc.log" 2>&1 || true
fi

####################
### peak calling and ChIP/Input ratio tracks
call_peaks_for_chip() {
    local sid=$1
    local ctrl=${control_by_id[$sid]}
    local chip_bam="$output_path/alignments/$sid/${sid}.filtered.bam"
    local ctrl_bam="$output_path/alignments/$ctrl/${ctrl}.filtered.bam"
    local peak_dir="$output_path/peaks/$sid"
    local paired=${paired_by_id[$sid]}
    local format="BAM"
    [[ "$paired" == "false" ]] || format="BAMPE"

    mkdir -p "$peak_dir"

    bamCompare \
        -b1 "$chip_bam" -b2 "$ctrl_bam" \
        -o "$output_path/tracks/${sid}.vs.${ctrl}.log2ratio.bw" \
        --operation log2 \
        --scaleFactorsMethod readCount \
        --pseudocount 1 1 \
        --binSize "$bin_size" \
        --numberOfProcessors "$cores_per_job" \
        > "$output_path/logs/${sid}.bamCompare.log" 2>&1

    local macs_args=(
        callpeak
        -t "$chip_bam"
        -c "$ctrl_bam"
        -f "$format"
        -g "$effective_genome_size"
        -n "$sid"
        --outdir "$peak_dir"
        --keep-dup all
    )

    if [[ "$broad_peaks" == "true" ]]; then
        macs_args+=(-q "$macs_qvalue" --broad --broad-cutoff "$broad_cutoff")
    else
        macs_args+=(-q "$macs_qvalue")
    fi

    macs3 "${macs_args[@]}" > "$output_path/logs/${sid}.macs3.log" 2>&1
}

log "Calling peaks with MACS3"
batch_pids=()
batch_labels=()
for ((u=0; u<${#chip_ids[@]}; u++)); do
    sid=${chip_ids[$u]}
    call_peaks_for_chip "$sid" > "$output_path/logs/${sid}.peak_pipeline.log" 2>&1 &
    batch_pids+=("$!")
    batch_labels+=("$sid")

    if (( ${#batch_pids[@]} == chunks )); then
        wait_for_batch batch_pids batch_labels
        batch_pids=()
        batch_labels=()
    fi
done
if (( ${#batch_pids[@]} > 0 )); then
    wait_for_batch batch_pids batch_labels
fi

####################
### build DiffBind sample sheet and run differential analysis
if [[ "$run_diff" == "true" ]]; then
    log "Preparing DiffBind sample sheet"
    diffbind_sheet="$output_path/diffbind/diffbind_samples.csv"
    echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,bamControl,Peaks,PeakCaller" > "$diffbind_sheet"

    compared_chip_count=0
    for sid in "${chip_ids[@]}"; do
        cond=${condition_by_id[$sid]}
        [[ "$cond" == "$test_condition" || "$cond" == "$reference_condition" ]] || continue

        ctrl=${control_by_id[$sid]}
        fac=${factor_by_id[$sid]}
        rep=${replicate_by_id[$sid]}
        chip_bam="$output_path/alignments/$sid/${sid}.filtered.bam"
        ctrl_bam="$output_path/alignments/$ctrl/${ctrl}.filtered.bam"

        if [[ "$broad_peaks" == "true" ]]; then
            peak_file="$output_path/peaks/$sid/${sid}_peaks.broadPeak"
            peak_caller="bed"
        else
            peak_file="$output_path/peaks/$sid/${sid}_peaks.narrowPeak"
            peak_caller="narrow"
        fi
        [[ -s "$peak_file" ]] || fail "MACS3 peak file is absent or empty for '$sid': $peak_file"

        printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "$sid" "Arabidopsis" "$fac" "$cond" "$rep" \
            "$chip_bam" "$ctrl_bam" "$peak_file" "$peak_caller" \
            >> "$diffbind_sheet"
        ((compared_chip_count+=1))
    done
    (( compared_chip_count >= 4 )) || fail "Too few ChIP samples were added to the DiffBind sheet."

    diffbind_r="$output_path/diffbind/run_diffbind.R"
    cat > "$diffbind_r" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
    stop("Expected 9 arguments: sheet outdir test reference minOverlap fdr lfc broad ncores")
}

sheet       <- normalizePath(args[1], mustWork = TRUE)
outdir      <- normalizePath(args[2], mustWork = TRUE)
test        <- args[3]
reference   <- args[4]
min_overlap <- as.integer(args[5])
fdr_cutoff  <- as.numeric(args[6])
lfc_cutoff  <- as.numeric(args[7])
broad       <- as.logical(args[8])
ncores      <- as.integer(args[9])

suppressPackageStartupMessages(library(DiffBind))
suppressPackageStartupMessages(library(BiocParallel))

if (.Platform$OS.type == "windows") {
    BiocParallel::register(BiocParallel::SnowParam(workers = ncores))
} else {
    BiocParallel::register(BiocParallel::MulticoreParam(workers = ncores))
}

if (is.na(min_overlap) || min_overlap < 1) stop("Invalid minOverlap")
if (is.na(fdr_cutoff) || fdr_cutoff <= 0 || fdr_cutoff > 1) stop("Invalid FDR threshold")
if (is.na(lfc_cutoff) || lfc_cutoff < 0) stop("Invalid LFC threshold")

options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = ncores)

message("Loading DiffBind sample sheet: ", sheet)
dba_obj <- dba(sampleSheet = sheet)

pdf(file.path(outdir, "01_peakset_correlation_heatmap.pdf"), width = 8, height = 8)
try(plot(dba_obj), silent = TRUE)
dev.off()

# Broad histone domains should not be re-centered to fixed-width summits.
if (broad) {
    dba_obj <- dba.count(
        dba_obj,
        minOverlap = min_overlap,
        summits = FALSE,
        bRemoveDuplicates = FALSE,
        bParallel = TRUE
    )
} else {
    dba_obj <- dba.count(
        dba_obj,
        minOverlap = min_overlap,
        summits = 250,
        bRemoveDuplicates = FALSE,
        bParallel = TRUE
    )
}

saveRDS(dba_obj, file.path(outdir, "diffbind_counted.rds"))

pdf(file.path(outdir, "02_count_correlation_heatmap.pdf"), width = 8, height = 8)
try(plot(dba_obj), silent = TRUE)
dev.off()

pdf(file.path(outdir, "03_PCA.pdf"), width = 8, height = 7)
try(dba.plotPCA(dba_obj, attributes = DBA_CONDITION), silent = TRUE)
dev.off()

message("Contrast: ", test, " vs ", reference)
dba_obj <- dba.contrast(
    dba_obj,
    contrast = c("Condition", test, reference),
    minMembers = 2
)

dba_obj <- dba.analyze(
    dba_obj,
    method = DBA_DESEQ2,
    bBlacklist = FALSE,
    bGreylist = FALSE,
    bParallel = TRUE
)
saveRDS(dba_obj, file.path(outdir, "diffbind_analyzed.rds"))

capture.output(
    dba.show(dba_obj, bContrasts = TRUE, th = fdr_cutoff),
    file = file.path(outdir, "contrast_summary.txt")
)

all_gr <- dba.report(
    dba_obj,
    contrast = 1,
    method = DBA_DESEQ2,
    th = 1,
    fold = 0,
    DataType = DBA_DATA_GRANGES
)
all_df <- as.data.frame(all_gr)
all_df <- all_df[order(all_df$FDR, -abs(all_df$Fold)), , drop = FALSE]
write.csv(all_df, file.path(outdir, "differential_peaks_all.csv"), row.names = FALSE)

sig_df <- all_df[
    !is.na(all_df$FDR) & all_df$FDR <= fdr_cutoff &
        !is.na(all_df$Fold) & abs(all_df$Fold) >= lfc_cutoff,
    , drop = FALSE
]
write.csv(sig_df, file.path(outdir, "differential_peaks_significant.csv"), row.names = FALSE)

write_bed <- function(df, filename) {
    if (nrow(df) == 0) {
        file.create(filename)
        return(invisible(NULL))
    }
    bed <- data.frame(
        chrom = as.character(df$seqnames),
        start = pmax(0, as.integer(df$start) - 1L),
        end = as.integer(df$end),
        name = paste0("DB_peak_", seq_len(nrow(df))),
        score = pmin(1000L, as.integer(round(-10 * log10(pmax(df$FDR, 1e-100))))),
        strand = ".",
        log2FoldChange = df$Fold,
        FDR = df$FDR,
        stringsAsFactors = FALSE
    )
    write.table(
        bed,
        filename,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
}
write_bed(sig_df, file.path(outdir, "differential_peaks_significant.bed"))

pdf(file.path(outdir, "04_MA_plot.pdf"), width = 8, height = 7)
try(dba.plotMA(dba_obj, contrast = 1, method = DBA_DESEQ2, th = fdr_cutoff), silent = TRUE)
dev.off()

pdf(file.path(outdir, "05_volcano_plot.pdf"), width = 8, height = 7)
try(dba.plotVolcano(dba_obj, contrast = 1, method = DBA_DESEQ2, th = fdr_cutoff), silent = TRUE)
dev.off()

summary_lines <- c(
    paste0("Contrast: ", test, " vs ", reference),
    paste0("Positive Fold: greater enrichment in ", test),
    paste0("Consensus minOverlap: ", min_overlap),
    paste0("FDR threshold: ", fdr_cutoff),
    paste0("Absolute log2 fold-change threshold: ", lfc_cutoff),
    paste0("All tested consensus sites: ", nrow(all_df)),
    paste0("Significant differential sites: ", nrow(sig_df)),
    paste0("Higher in ", test, ": ", sum(sig_df$Fold > 0, na.rm = TRUE)),
    paste0("Higher in ", reference, ": ", sum(sig_df$Fold < 0, na.rm = TRUE))
)
writeLines(summary_lines, file.path(outdir, "differential_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "R_sessionInfo.txt"))
RSCRIPT

    log "Running DiffBind/DESeq2"
    Rscript "$diffbind_r" \
        "$diffbind_sheet" \
        "$output_path/diffbind" \
        "$test_condition" \
        "$reference_condition" \
        "$min_overlap" \
        "$fdr_threshold" \
        "$lfc_threshold" \
        "$broad_peaks" \
        "$n_cores" \
        > "$output_path/logs/diffbind.log" 2>&1
fi

####################
### final summary
{
    echo "Pipeline completed: $(date -Is)"
    echo "Output directory: $output_path"
    echo "Effective genome size: $effective_genome_size"
    echo "Peak mode: $([[ "$broad_peaks" == "true" ]] && echo broad || echo narrow)"
    echo "ChIP samples: ${chip_ids[*]}"
    if [[ "$run_diff" == "true" ]]; then
        echo "Contrast: $test_condition vs $reference_condition"
        echo "DiffBind results: $output_path/diffbind"
    fi
} > "$output_path/run_summary.txt"

rm -rf "$output_path/tmp"
log "Pipeline completed successfully"
echo "**  results: $output_path"
