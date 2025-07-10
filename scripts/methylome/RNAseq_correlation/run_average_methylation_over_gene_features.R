sse_h_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_240325/Methylome.At/samples_table/samples_table_sseHigh.txt"
sse_l_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_240325/Methylome.At/samples_table/samples_table_sseLow.txt"
sse_h_l_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_240325/Methylome.At/samples_table/samples_table_SSE_high_low.txt"
dcgs_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_240325/Methylome.At/samples_table/samples_table_dCGS.txt"

output_path <- "/home/yoyerush/yo/methylome_pipeline/"

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/methylome/RNAseq_correlation/sample_ave_methylation_levels_from_gene_features.R")
average_methylation_over_gene_features(sse_h_path, output_path)
average_methylation_over_gene_features(sse_l_path, output_path)
average_methylation_over_gene_features(sse_h_l_path, output_path)
average_methylation_over_gene_features(dcgs_path, output_path)
#