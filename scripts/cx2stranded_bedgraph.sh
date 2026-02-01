#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   cx2stranded_bedgraph.sh input.CX_report.txt.gz output.bedGraph
#   cx2stranded_bedgraph.sh input.CX_report.txt.gz            # output auto

in="${1:-}"
out="${2:-}"

if [[ -z "$in" || ! -f "$in" ]]; then
    echo "Usage: $(basename "$0") input.CX_report.txt.gz [output.bedGraph]"
    exit 1
fi

# Default output name if not provided
if [[ -z "$out" ]]; then
    out="${in%.gz}.stranded.bedGraph"
    out="${out%.txt}.stranded.bedGraph"
fi

# Use gunzip -c if .gz, else cat
if [[ "$in" == *.gz ]]; then
    reader="gunzip -c"
else
    reader="cat"
fi

$reader "$in" \
| awk 'BEGIN{OFS="\t"}
    {
        m=$4; u=$5; cov=m+u;
        if(cov==0) next;
        perc=100*m/cov;
        start=$2-1; end=$2;
        if($3=="+") print $1,start,end, perc;
        else if($3=="-") print $1,start,end,-perc;
    }' \
| sort -k1,1 -k2,2n > "$out"

echo "Wrote: $out"