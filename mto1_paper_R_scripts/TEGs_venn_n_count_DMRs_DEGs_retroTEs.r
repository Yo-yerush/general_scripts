library(dplyr)
library(VennDiagram)

############################################

DMR_file.0 <- rbind(
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/TEG_CG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/TEG_CHG_genom_annotations.csv"),
    read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/TEG_CHH_genom_annotations.csv")
)

DMR_file <- DMR_file.0 %>%
    select(gene_id, regionType, context) %>%
    distinct(gene_id, .keep_all = T)

############################################

RNA_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    filter(gene_model_type == "transposable_element_gene") %>%
    dplyr::rename(gene_id = "locus_tag") %>%
    filter(padj < 0.05) %>%
    select(gene_id, log2FoldChange) %>%
    distinct(gene_id, .keep_all = T)

############################################

# use TEGs with 'Derives_from' column
gene_2_TE_ids <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_description_file.csv.gz") %>%
    filter(!is.na(Derives_from)) %>%
    distinct(gene_id, Derives_from)

# filter to retro-TEs
TE_file <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt",
    sep = "\t"
) %>%
    mutate(seqnames = NA) %>% # Add a new column with NA values
    dplyr::rename(Derives_from = Transposon_Name) %>%
    dplyr::select(Derives_from, Transposon_Super_Family) %>%
    filter(grepl("Gypsy|Copia|LINE", Transposon_Super_Family)) %>%
    merge(., gene_2_TE_ids, by = "Derives_from")

############################################

retro_DMRs <- merge(DMR_file, TE_file, by = "gene_id")
retro_DEGs <- merge(RNA_file, TE_file, by = "gene_id")

############################################
gene_sets <- list(
    DMRs = retro_DMRs$gene_id,
    DEGs = retro_DEGs$gene_id
)

# venn_colors <- c("#a05b9c", "#71c071")
venn_colors <- c("#928e92", "#d69641")
category.position <- c(0, 0)
resolution <- 300
cex <- 0.75

venn.diagram(
    x = gene_sets,
    category.names = c(
        "\nunique TEGs                    \noverlap DMRs                   ",
        "\n\n\n\n\n\n\n\nDE-TEGs"
    ),
    filename = paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/DMRs_over_de-retro-TEGs_VennDiagram.png"),
    disable.logging = T,
    output = T,
    imagetype = "png",
    height = 480,
    width = 480,
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
    cat.pos = category.position,
    cat.fontface = 2,
    cat.fontfamily = "serif"
    #    cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
    #    col=venn_colors,
    #    rotation = 1
)

message(
    "\n\n",
    "hyper-DMRs: ",
    (nrow(filter(retro_DMRs, regionType == "gain")) / nrow(retro_DMRs)) * 100,
    "\n",
    "downregulated DEGs: ",
    (nrow(filter(retro_DEGs, log2FoldChange < 0)) / nrow(retro_DEGs)) * 100,
    "\n\n",
    "overlap DMRs (unique TEGs): ",
    nrow(retro_DMRs),
    "\n",
    "overlap DMRs (not unique TEGs): ",
    DMR_file.0 %>%
        merge(., TE_file, by = "gene_id") %>%
        mutate(tmp = paste(seqnames, start, end, sep = "_")) %>%
        distinct(tmp) %>%
        nrow(),
    "\n\n",
    "DE-TEGs: ",
    nrow(retro_DEGs),
    "\n",
    "overlap DE-TEGs: ",
        length(intersect(gene_sets[[1]], gene_sets[[2]])),
    "\n"
)
