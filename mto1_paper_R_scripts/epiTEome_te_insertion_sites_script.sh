#!/bin/bash

#-------------------------------------------------------------
#
### first run:
# ./run_bismark_yo.sh -s ./mto1_samples_bismark.txt -g /home/yoyerush/yo/methylome_pipeline/Bismark/res_unmapped_040825/TAIR10_chr_all.fas.gz -o bismark_results -n 32 -m 16G --um
#
#-------------------------------------------------------------
#
### run to create retro-TEs list:
# awk '$9 ~ /sF=LTR\/Copia|sF=LTR\/Gypsy|sF=LINE\/L1/ {
#     split($9, arr, ";");
#     for (i in arr) {
#         if (arr[i] ~ /ID=/) {
#             split(arr[i], id, "=");
#             print id[2];
#         }
#     }
# }' tair10TEs.gff3 | sort | uniq > retrotransposons_list.txt
#
### run to create helitron - "ATREP9" family list:
# awk '$9 ~ /fam=ATREP9/ {
#     split($9, arr, ";");
#     for (i in arr) {
#         if (arr[i] ~ /ID=/) {
#             split(arr[i], id, "=");
#             print id[2];
#         }
#     }
# }' tair10TEs.gff3 | sort | uniq > helitron_atrep9_list.txt
#
### run to create copia superfamily list:
# awk '$9 ~ /sF=LTR\/Copia/ {
#     split($9, arr, ";");
#     for (i in arr) {
#         if (arr[i] ~ /ID=/) {
#             split(arr[i], id, "=");
#             print id[2];
#         }
#     }
# }' tair10TEs.gff3 | sort | uniq > copia_list.txt
#
### run to create gypsy superfamily list:
# awk '$9 ~ /sF=LTR\/Gypsy/ {
#     split($9, arr, ";");
#     for (i in arr) {
#         if (arr[i] ~ /ID=/) {
#             split(arr[i], id, "=");
#             print id[2];
#         }
#     }
# }' tair10TEs.gff3 | sort | uniq > gypsy_list.txt
#
### run to create LIN!/L1 superfamily list:
# awk '$9 ~ /sF=LINE\/L1/ {
#     split($9, arr, ";");
#     for (i in arr) {
#         if (arr[i] ~ /ID=/) {
#             split(arr[i], id, "=");
#             print id[2];
#         }
#     }
# }' tair10TEs.gff3 | sort | uniq > line_l1_list.txt
#
#-------------------------------------------------------------

# epiTEome scripts dir - also the output files path
cd /home/yoyerush/yo/methylome_pipeline/transposition_activity/

samples_name=("mto1_1" "mto1_2" "mto1_3" "wt_1" "wt_2")

te_ids_list=copia_list.txt

bismark_results=/home/yoyerush/yo/methylome_pipeline/Bismark/res_unmapped_040825/bismark_results

#-------------------------------------------------------------

# index the genome
perl ./idxEpiTEome.pl -l 150 -gff tair10TEs.gff3 -t $te_ids_list -fasta TAIR10_chr_all.fna

mkdir -p unmapped_reads

# Loop over samples
for sample in "${samples_name[@]}"; do
    # concatenate paired unmapped reads
    zcat $bismark_results/"$sample"/"$sample"_unmapped_reads_1.fq.gz $bismark_results/"$sample"/"$sample"_unmapped_reads_2.fq.gz > unmapped_reads/"$sample"_unmapped_reads.fq
    
    # run epiTEome
    perl ./epiTEome.pl -gff tair10TEs.gff3 -t $te_ids_list -ref TAIR10_chr_all.epiTEome.masked.fasta -un unmapped_reads/"$sample"_unmapped_reads.fq -p 20
    
done

# gzip unmapped_reads/*_unmapped_reads.fq
