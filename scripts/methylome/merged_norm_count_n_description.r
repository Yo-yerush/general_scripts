library(dplyr)

ann_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_description_file.csv.gz") %>%
    select(gene_id, Symbol, Computational_description)

RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/norm.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = X)

#######################

merged = merge.data.frame(RNAseq, ann_file, by = "gene_id", all.x = T)
merged <- merged %>% mutate(across(where(is.numeric), round, digits = 1)) # round

wt_columns = grep("met20|met21|met22", names(merged))
mto1_columns = grep("met14|met15|met16", names(merged))

names(merged)[wt_columns] <- paste0("wt_", 1:3)
names(merged)[mto1_columns] <- paste0("mto1_", 1:3)

merged$mto1_sum = rowSums(merged[mto1_columns], na.rm = TRUE)
merged$wt_sum = rowSums(merged[wt_columns], na.rm = TRUE)

merged <- arrange(merged, wt_sum) %>%
    dplyr::relocate(Symbol, .after = gene_id) %>%
    dplyr::relocate(Computational_description, .after = wt_sum)

merged$Symbol[is.na(merged$Symbol)] = ""

openxlsx::write.xlsx(
    merged,
    "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/norm_count_n_description.xlsx"
)
