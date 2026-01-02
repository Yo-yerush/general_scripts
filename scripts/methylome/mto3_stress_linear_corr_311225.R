var1_name <- "wt"
var2_name <- "mto3"
var1_rnaseq_names <- c("met20", "met22")
var2_rnaseq_names <- c("met17", "met18", "met19")

average_meth_results_directory <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/average.meth.genes.levels/mto3_vs_wt/"
DMRs_results_directory <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto3_vs_wt/genome_annotation/"
RNAseq_results_directory <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto3_vs_wt/"
main_output_directory <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/Linear_correlation/"

#-------------------------------------------------------------------#

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/multiplot_ggplot2.R")
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/yo_theme_base_ggplot2.r")
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/linear_plot_over_gene_features_DMRs_DEGs.r")

#-------------------------------------------------------------------#

offspring_fun <- function(go_id, xx = as.list(GO.db::GOBPOFFSPRING)) { # 'GOBPCHILDREN' for child terms

    child_terms_0 <- as.character(xx[[go_id]])
    child_terms <- child_terms_0

    for (i in 1:length(child_terms_0)) {
        child_terms <- c(child_terms, as.character(xx[[child_terms[i]]]))
    }

    child_terms <- child_terms[!is.na(child_terms)] %>% unique()

    out <- AnnotationDbi::select(
        x = org.At.tair.db,
        keys = child_terms,
        columns = c("GO", "TAIR"),
        keytype = "GO"
    ) %>%
        filter(!is.na(TAIR)) %>%
        distinct(TAIR)

    c(out)$TAIR
}

library(dplyr)
library(org.At.tair.db)

child_terms_stress <- offspring_fun("GO:0006950") # response to stress
child_terms_biotic <- offspring_fun("GO:0009607") # response to biotic stimulus
child_terms_abiotic <- offspring_fun("GO:0009628") # response to abiotic stimulus

#-------------------------------------------------------------------#

run_list <- list(
    stress = list(child_terms_stress, "#2e720e"),
    biotic = list(child_terms_biotic, "#a12622"),
    abiotic = list(child_terms_abiotic, "#6e22a1")
)

for (i_run in c("stress", "biotic", "abiotic")) {
    linear_plot_meth_rna(
        var1_name = var1_name,
        var2_name = var2_name,
        var1_rnaseq_names = var1_rnaseq_names,
        var2_rnaseq_names = var2_rnaseq_names,
        average_meth_results_directory = average_meth_results_directory,
        DMRs_results_directory = DMRs_results_directory,
        RNAseq_results_directory = RNAseq_results_directory,
        genes_2_keep = run_list[[i_run]][[1]],
        genes_group_name = i_run,
        main_output_directory = main_output_directory,
        additional_plots = TRUE,
        pValues_table = TRUE,
        var2_col = run_list[[i_run]][[2]]
    )
}



# var1_name = var1_name
# var2_name = var2_name
# var1_rnaseq_names = var1_rnaseq_names
# var2_rnaseq_names = var2_rnaseq_names
# average_meth_results_directory = average_meth_results_directory
# DMRs_results_directory = DMRs_results_directory
# RNAseq_results_directory = RNAseq_results_directory
# genes_2_keep = run_list[[1]][[1]]
# main_output_directory = main_output_directory
# additional_plots = TRUE
# pValues_table = TRUE
# var2_col = run_list[[1]][[2]]