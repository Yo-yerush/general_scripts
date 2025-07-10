
sample_file_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_240325/Methylome.At/samples_table/samples_table_mto1.txt"
output_path <- "/home/yoyerush/yo/methylome_pipeline/"
n.cores <- 20

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/methylome/RNAseq_correlation/sample_ave_methylation_levels_from_gene_features.R")
average_methylation_over_gene_features(sample_file_path, output_path, n.cores)