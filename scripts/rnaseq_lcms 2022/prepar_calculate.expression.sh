#!/bin/bash

cd /home/yoyerush/yo/all/index
rsem-prepare-reference --bowtie2 /home/yoyerush/yo/Reference_genome_ASM765513v2/ncbi_dataset/data/GCF_007655135.1/rna.fna indx

for i in {1..30}
do
mkdir /home/yoyerush/yo/all/peel/expression/P"$i"
cd /home/yoyerush/yo/all/peel/expression/P"$i"
rsem-calculate-expression -p 10 --bowtie2 /home/yoyerush/yo/all/peel/trimming/P"$i"_trim.fastq /home/yoyerush/yo/all/index/indx PEEL"$i" --append-names
done
