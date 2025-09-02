library(dplyr)
library(VennDiagram)
library(org.At.tair.db)
library(KEGGREST)

path_2_save <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/שיט/וועדה דוקטורט טכניון_2025/"

down_mto1 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all_genes_results_mto1_vs_wt.csv") %>%
    filter(padj < 0.05) %>%
    filter(log2FoldChange < 0) %>%
    distinct(gene_id)

down_mto3 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto3_vs_wt/all_genes_results_mto3_vs_wt.csv") %>%
    filter(padj < 0.05) %>%
    filter(log2FoldChange < 0) %>%
    distinct(gene_id)

############################################

up_mto1 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all_genes_results_mto1_vs_wt.csv") %>%
    filter(padj < 0.05) %>%
    filter(log2FoldChange > 0) %>%
    distinct(gene_id)

up_mto3 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto3_vs_wt/all_genes_results_mto3_vs_wt.csv") %>%
    filter(padj < 0.05) %>%
    filter(log2FoldChange > 0) %>%
    distinct(gene_id)

############################################

venn_from_map <- function(map_id, down = T) {
    met_genes <- data.frame(gene_id = as.list(org.At.tairPATH2TAIR)[[map_id]])

    if (down) {
        gene_sets <- list(
            mto1 = down_mto1$gene_id,
            mto3 = down_mto3$gene_id,
            Met_pathway = met_genes$gene_id
        )
    } else {
        gene_sets <- list(
            mto1 = up_mto1$gene_id,
            mto3 = up_mto3$gene_id,
            Met_pathway = met_genes$gene_id
        )
    }

    # venn_colors <- c("#a05b9c", "#71c071")
    venn_colors <- c("#928e92", "#d69641", "#71c071")
    category.position <- c(0, 0)
    resolution <- 300
    cex <- 0.6

    venn.diagram(
        x = gene_sets,
        category.names = c("mto1", "mto3", paste0("ath", map_id)),
        filename = paste0(path_2_save, "venn_mtos_", ifelse(down, "down", "up"), "_met_pathway_ath", map_id, ".png"),
        disable.logging = T,
        output = T,
        imagetype = "png",
        height = 380,
        width = 380,
        resolution = resolution,
        compression = "lzw",
        lwd = 1,
        fill = venn_colors[1:length(gene_sets)],
        alpha = rep(0.45, length(gene_sets)),
        col = rep("white", length(gene_sets)),
        cex = cex,
        fontfamily = "serif",
        cat.cex = cex,
        cat.default.pos = "outer",
        # cat.pos = category.position,
        cat.fontface = 2,
        cat.fontfamily = "serif"
        #    cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
        #    col=venn_colors,
        #    rotation = 1
    )
}

for (t_or_F in c(T,F)) {
    venn_from_map("00270", t_or_F) # Cysteine and methionine metabolism
    venn_from_map("00920", t_or_F) # Sulfur metabolism
    venn_from_map("00966", t_or_F) # Glucosinolate biosynthesis
    venn_from_map("00500", t_or_F) # Starch and sucrose metabolism
}
