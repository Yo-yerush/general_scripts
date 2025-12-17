library(dplyr)


des_file <- read.csv("C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/RA_costume_annotations_files/Methylome.At_description_file.csv") %>%
select(gene_id, Symbol, Short_description)

load_df <- function(x) {
    xx <- readxl::read_excel(paste0("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/genes_group_results/",x,"_groups.xlsx"),
    sheet = "seed_specific_genes") %>%
    filter(RNA_pvalue < 0.05, RNA_log2FC > 0, r_value >= 0.7) %>%
    select(gene_id,RNA_log2FC,RNA_pvalue)
    names(xx)[-1] <- paste(names(xx)[-1], x, sep = "_")
    names(xx) <- gsub("RNA_", "", names(xx))
    xx
}

mto1 <- load_df("mto1_vs_wt")
mto3 <- load_df("mto3_vs_wt")
hm <- load_df("SSE_high_vs_EV")

merged_df <- merge(mto1, mto3) %>% merge(., hm) %>% merge(., des_file)
writexl::write_xlsx(merged_df, "C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_seed_specific_upregulation.xlsx")
