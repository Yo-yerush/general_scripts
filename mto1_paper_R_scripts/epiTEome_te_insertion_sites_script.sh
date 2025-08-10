#!/bin/bash
#-------------------------------------------------------------
# EpiTEome tool
# https://github.com/jdaron/epiTEome
# https://doi.org/10.1186/s13059-017-1232-0
#-------------------------------------------------------------
#
# Usage: ./epiTEome_te_insertion_sites_script.sh [--dont_indx]
#   --dont_indx          Skip genome indexing step
#   --dont_concatenate   Skip concatenation of unmapped reads
#   --run_test           Use test data
#
#-------------------------------------------------------------
#
# # first get unmapped reads using:
# ./run_bismark_yo.sh -s ./mto1_samples_bismark.txt -g /home/yoyerush/yo/methylome_pipeline/Bismark/res_unmapped_040825/TAIR10_chr_all.fas.gz -o bismark_results -n 32 -m 16G --um
#
#-------------------------------------------------------------
#
# # if there is error writing: << Unknown command 'filter' >>
# # then modify the 'bamutils' script. example in 'install_rpiTEome_env.sh' script
#
#-------------------------------------------------------------

# epiTEome scripts dir - also the output files path
cd /home/yoyerush/yo/methylome_pipeline/transposition_activity
samples_name=("mto1_1" "mto1_2" "mto1_3" "wt_1" "wt_2")
te_super_family=copia
te_ids_list=te_lists/copia_list.txt
bismark_results=/home/yoyerush/yo/methylome_pipeline/Bismark/res_unmapped_040825/bismark_results
genome_file=genome_indx/TAIR10_chr_all.epiTEome.masked.fasta
gff_file=genome_indx/tair10TEs.gff3
n_cores=30

#-------------------------------------------------------------

# '--dont_concatenate' arguments
dont_concatenate=0
for arg in "$@"; do
    if [[ "$arg" == "--dont_concatenate" ]]; then
        dont_concatenate=1
        break
    fi
done

# '--dont_indx' arguments
dont_indx=0
for arg in "$@"; do
    if [[ "$arg" == "--dont_indx" ]]; then
        dont_indx=1
        break
    fi
done

# '--run_test' arguments
run_test=0
for arg in "$@"; do
    if [[ "$arg" == "--run_test" ]]; then
        run_test=1
        break
    fi
done

#-------------------------------------------------------------

mkdir -p results

if [[ $run_test -eq 0 ]]; then
    
    # index the genome
    if [[ $dont_indx -eq 0 ]]; then
        perl epiteome_scripts/idxEpiTEome.pl -l 150 -gff genome_indx/tair10TEs.gff3 -t $te_ids_list -fasta genome_indx/TAIR10_chr_all.fna
    fi
    
    mkdir -p unmapped_reads
    
    # Loop over samples
    for sample in "${samples_name[@]}"; do
        
        # concatenate paired unmapped reads
        if [[ $dont_concatenate -eq 0 ]]; then
            zcat $bismark_results/"$sample"/"$sample"_unmapped_reads_1.fq.gz $bismark_results/"$sample"/"$sample"_unmapped_reads_2.fq.gz > unmapped_reads/"$sample"_unmapped_reads.fq
        fi
        
        mkdir -p results/"$sample"
        cp unmapped_reads/"$sample"_unmapped_reads.fq results/"$sample"
        
        # run epiTEome
        perl epiteome_scripts/epiTEome_yo_2.pl -gff $gff_file -t $te_ids_list -ref $genome_file -un results/"$sample"/"$sample"_unmapped_reads.fq -p $n_cores
        
        
        mkdir -p results_"$te_super_family"
        ### fix this shit ###
        ### mv unmapped_reads/"$sample"unmapped.newInsertionSite.tab, unmapped.newInsertionSite.sam, unmapped.met.meta.tab and unmapped.met.row.tab
    done
    
else
    mkdir results/test_data_res
    mkdir results/test_data_res/genome_indx_test
    cp test_data/unmapped.fastq.bz2 results/test_data_res
    cp test_data/Chr2.fasta results/test_data_res/genome_indx_test
    perl epiteome_scripts/idxEpiTEome.pl -l 150 -gff genome_indx/tair10TEs.gff3 -t test_data/teid.lst -fasta results/test_data_res/genome_indx_test/Chr2.fasta
    perl epiteome_scripts/epiTEome_yo_2.pl -gff $gff_file -t test_data/teid.lst -ref results/test_data_res/genome_indx_test/Chr2.fasta -un results/test_data_res/unmapped.fastq.bz2 -p 1
fi