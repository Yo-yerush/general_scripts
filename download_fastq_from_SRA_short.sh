#!/bin/bash

# Usage: 
#   download_fastq_from_SRA.sh "SRR534177 SRR534193 SRR534209" "/path/to/output"

SRR_ARRAY=($1) # space-separated SRR IDs
OUTPUT_DIR=$2

########################

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

for SRR_ID in "${SRR_ARRAY[@]}"; do
    echo "Downloading ${SRR_ID}..."
    prefetch "$SRR_ID"
    
    # convert the .sra file to fastq.gz
    echo "Converting ${SRR_ID} to FASTQ..."
    fasterq-dump "$SRR_ID" --split-files --threads 10 # --outdir "$OUTPUT_DIR"
    gzip "$SRR_ID".fastq
    
    # remove the .sra file
    echo "Removing ${SRR_ID}.sra file..."
    rm -r "${SRR_ID}"
    echo ""
done

########################

echo "completed. check '$OUTPUT_DIR' directory"