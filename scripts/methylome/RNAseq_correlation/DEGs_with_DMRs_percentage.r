
library(dplyr)

final_df = data.frame(DMRs.position = c("Gene.body", "Promoter"), CG = NA, non.CG = NA, All.contexts = NA)

DEGs <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag) %>%
    filter(padj < 0.05) %>%
    dplyr::select(gene_id)

percentage_fun <- function(col.name, ann) {
    x_cg <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/", ann, "_CG_genom_annotations.csv"))

    x_chg <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/", ann, "_CHG_genom_annotations.csv"))

    x_chh <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/", ann, "_CHH_genom_annotations.csv"))

    if (col.name == "CG") {
        x = x_cg

    } else if(col.name == "non.CG") {
        x = rbind(x_chg, x_chh)

    } else if(col.name == "All.contexts") {
        x = rbind(x_cg, x_chg, x_chh)
    }
    
    x <- x %>%
        distinct(gene_id) %>%
        merge.data.frame(DEGs) %>%
        nrow() / nrow(DEGs)

    xx <- paste0(round(x * 100, 2), "%")

    return(xx)
}

for (i_column in c("CG", "non.CG", "All.contexts")) {
   final_df[1, i_column] <- percentage_fun(i_column, "Genes")
   final_df[2, i_column] <- percentage_fun(i_column, "Promoters")
}
