#!/bin/bash
#-------------------------------------------------------------
# epiTEome tool
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
# # there is error writing: << Unknown command 'filter' >>
# # using midofied epiTEome.pl script - use samtools=1.16 instead of the 'ngsutils' package
# # also can modify the 'bamutils' script. example in 'install_3piTEome_env.sh' script
#
#-------------------------------------------------------------

# epiTEome scripts dir - also the output files path
main_dir=/home/yoyerush/yo/methylome_pipeline/transposition_epiTEome
samples_name=("mto1_1" "mto1_2" "mto1_3" "wt_1" "wt_2")
te_super_family=copia
te_ids_list=copia_list.txt
bismark_results=/home/yoyerush/yo/methylome_pipeline/Bismark/res_unmapped_040825/bismark_results
genome_file=TAIR10_chr_all.fna
gff_file=tair10TEs.gff3
n_cores=30

# genome indexed file with - '.epiTEome.masked.fasta' suffix
genome_indx_file=$(basename "$genome_file" | sed 's/\.[^.]*$/.epiTEome.masked.fasta/')

cd $main_dir

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

mkdir -p "$main_dir"/results

if [[ $run_test -eq 0 ]]; then
    
    cd $main_dir

    # index the genome
    if [[ $dont_indx -eq 0 ]]; then
        perl epiteome_scripts/idxEpiTEome.pl -l 150 -gff genome_indx/$gff_file -t te_lists/$te_ids_list -fasta genome_indx/$genome_file
    fi
    
    mkdir -p "$main_dir"/unmapped_reads
    
    # Loop over samples
    for sample in "${samples_name[@]}"; do
        
        # concatenate paired unmapped reads
        if [[ $dont_concatenate -eq 0 ]]; then
            zcat $bismark_results/"$sample"/"$sample"_unmapped_reads_1.fq.gz $bismark_results/"$sample"/"$sample"_unmapped_reads_2.fq.gz > "$main_dir"/unmapped_reads/"$sample"_unmapped_reads.fq
        fi
        
        mkdir -p "$main_dir"/results/"$sample"
        cp "$main_dir"/unmapped_reads/"$sample"_unmapped_reads.fq "$main_dir"/results/"$sample"
        cd "$main_dir"/results/"$sample"

        # run epiTEome
        perl ../../epiteome_scripts/epiTEome_Yo_edit.pl -gff ../../genome_indx/"$gff_file" -t ../../te_lists/$te_ids_list -ref ../../genome_indx/$genome_indx_file -un "$sample"_unmapped_reads.fq -p $n_cores
        
        cd $main_dir

        # # mkdir -p results_"$te_super_family"
        # mkdir -p "$main_dir"/results_"$te_super_family"/"$sample"
        # mv "$main_dir"/results/"$sample"/"$sample"_unmapped_reads.newInsertionSite.tab "$main_dir"/results_"$te_super_family"/"$sample"/
        # mv "$main_dir"/results/"$sample"/"$sample"_unmapped_reads.newInsertionSite.sam "$main_dir"/results_"$te_super_family"/"$sample"/
        # mv "$main_dir"/results/"$sample"/"$sample"_unmapped_reads.met.meta.tab "$main_dir"/results_"$te_super_family"/"$sample"/
        # mv "$main_dir"/results/"$sample"/"$sample"_unmapped_reads.met.row.tab "$main_dir"/results_"$te_super_family"/"$sample"/

        echo ""
        echo "***   Completed processing $sample    ***"
        echo ""
    done
    
else
    mkdir "$main_dir"/results/test_data_res
    mkdir "$main_dir"/results/test_data_res/genome_indx_test
    cp "$main_dir"/test_data/unmapped.fastq.bz2 "$main_dir"/results/test_data_res
    cp "$main_dir"/test_data/Chr2.fasta "$main_dir"/results/test_data_res/genome_indx_test

    cd "$main_dir"/results/test_data_res

    perl ../../epiteome_scripts/idxEpiTEome.pl -l 150 -gff ../../genome_indx/tair10TEs.gff3 -t ../../test_data/teid.lst -fasta genome_indx_test/Chr2.fasta

    cd "$main_dir"/results/test_data_res

    perl ../../epiteome_scripts/epiTEome_Yo_edit.pl -gff ../../genome_indx/tair10TEs.gff3 -t ../../test_data/teid.lst -ref genome_indx_test/Chr2.epiTEome.masked.fasta -un unmapped.fastq.bz2 -p 1

    cd $main_dir
fi