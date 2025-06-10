#!/bin/bash

usage() {
  echo "Usage: $0 [-t n_threads] [-o output_dir] SRR_ID [SRR_ID ...]"
  exit 1
}

# Set default values
n_threads=10
output_dir="./"

# Process optional arguments with getopts
while getopts "t:o:h" opt; do
  case "$opt" in
    t) n_threads="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

shift $((OPTIND - 1))

# Check if at least one SRR ID is provided
if [ "$#" -lt 1 ]; then
  echo ""
  echo "Error: Provide at least one SRR ID."
  echo ""
  usage
  echo ""
fi

# All remaining arguments are SRR IDs.
SRR_ARRAY=("$@")

echo "Using output directory: $output_dir"
echo "Using $n_threads thread(s) for fasterq-dump."
mkdir -p "$output_dir"
cd "$output_dir" || { echo "Could not change directory to $output_dir"; exit 1; }

#################################################

for SRR_ID in "${SRR_ARRAY[@]}"; do
    echo "Downloading ${SRR_ID}..."
    prefetch "$SRR_ID"
    
    # Convert .sra to FASTQ (using the specified number of threads)
    echo "Converting ${SRR_ID} to fastq..."
    fasterq-dump "$SRR_ID" --split-files --threads "$n_threads"
    gzip "${SRR_ID}"*.fastq
    
    # Remove the .sra directory
    echo "Removing ${SRR_ID}.sra file or directory..."
    rm -r "${SRR_ID}"
    echo ""
done

echo "Completed. Check '$output_dir' directory."
