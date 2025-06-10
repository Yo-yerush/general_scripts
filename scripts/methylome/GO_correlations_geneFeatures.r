library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(topGO)

##################
just_BP = FALSE
#################



RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag, DEG_log2FC = log2FoldChange) %>%
    dplyr::select(gene_id, DEG_log2FC, padj, pValue)

feature_file_fun <- function(context) {
    feature_file <- data.frame()
    for (ann in c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs", "TEG")) {
        # DMR results file
        ann_DMRs <- read.csv(paste0(
            "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/",
            context, "/", ann, "_", context, "_genom_annotations.csv"
        )) %>%
            dplyr::select(
                gene_id, log2FC, context, type, Symbol, Computational_description
            ) %>%
            dplyr::rename(DMR_log2FC = log2FC)

        # correlation results file
        ann_corr <- read.csv(paste0(
            "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/Gene_feature/mto1/", context, "/", ann, ".corr.", context, ".mto1.csv"
        )) %>%
            dplyr::select(-padj) %>%
            dplyr::rename(cor_pValue = pval)

        ann_merged <- merge.data.frame(ann_corr, ann_DMRs, by = "gene_id", all.y = T)

        feature_file <- rbind(feature_file, ann_merged)
    }
    return(feature_file)
}

################
remove_dup_DMR <- function(y) {
    y <- as.character(unique(unlist(strsplit(y, ","))))
    paste(y, collapse = ",")
}
################

################
feture_res_table <- rbind(
    feature_file_fun("CG"),
    feature_file_fun("CHG"),
    feature_file_fun("CHH")
) %>%
    mutate(
        DMR_log2FC = round(DMR_log2FC, 2),
        CG_DMRs = ifelse(grepl("CG", context), DMR_log2FC, NA),
        CHG_DMRs = ifelse(grepl("CHG", context), DMR_log2FC, NA),
        CHH_DMRs = ifelse(grepl("CHH", context), DMR_log2FC, NA)
    ) %>%
    merge.data.frame(RNAseq, ., by = "gene_id", all.y = TRUE) %>%
    mutate(tmp = paste(gene_id, type, sep = "_")) %>%
    group_by(tmp) %>%
    summarise(
        across(contains("CG_DMRs") | contains("CHG_DMRs") | contains("CHH_DMRs"), ~ remove_dup_DMR(paste(., collapse = ","))), # apply remove_dup_DMR function for DMR columns
        across(!contains("CG_DMRs") | contains("CHG_DMRs") | contains("CHH_DMRs"), dplyr::first) # for other columns
    ) %>%
    as.data.frame() %>%
    mutate(across(contains("_DMRs"), ~ gsub("NA", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",,", ",", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub("^,", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",$", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",", ", ", .))) %>%
    dplyr::relocate(CG_DMRs, CHG_DMRs, CHH_DMRs, .before = context) %>%
    dplyr::relocate(type, cor, cor_pValue, .after = gene_id) %>%
    dplyr::select(-context, -DMR_log2FC, -tmp) %>%
    # filter(pValue < 0.05) %>%
    arrange(pValue) %>%
    arrange(type) %>%
    mutate(across(contains("padj") | contains("pValue"), ~ gsub(" NA", NA, .)),
        pValue = as.numeric(formatC(.$pValue, format = "e", digits = 2)),
        padj = as.numeric(formatC(.$padj, format = "e", digits = 2)),
        DEG_log2FC = round(DEG_log2FC, 3),
        cor = round(cor, 3),
        cor_pValue = as.numeric(formatC(.$cor_pValue, format = "e", digits = 2))
    ) %>%
    ### add to this script:
    dplyr::select(gene_id, pValue, cor_pValue, DEG_log2FC) %>%
    distinct(gene_id, .keep_all = T)


feture_res_table$pValue[is.na(feture_res_table$pValue)] <- 0.999
feture_res_table$pValue[feture_res_table$pValue == 0] <- 1e-300



direction_fun <- function(x, gainORloss) {
    if (gainORloss == "gain") {
        geneList <- ifelse(x$pValue < 0.05 & x$cor_pValue < 0.05 & x$DEG_log2FC > 0, 1, 0) # filter by pValue, cor_pValue and up/down regulated
        names(geneList) <- x$gene_id
    } else if (gainORloss == "loss") {
        geneList <- ifelse(x$pValue < 0.05 & x$cor_pValue < 0.05 & x$DEG_log2FC < 0, 1, 0) # filter by pValue, cor_pValue and up/down regulated
        names(geneList) <- x$gene_id
    }

    res_list <- list()
    for (GO_type_loop in c("BP", "MF", "CC")) {
        myGOdata <- new("topGOdata",
            ontology = GO_type_loop,
            allGenes = geneList,
            geneSelectionFun = function(x) (x == 1),
            # description = "Test",
            annot = annFUN.org,
            # nodeSize = 5,
            mapping = "org.At.tair.db"
        )

        sg <- sigGenes(myGOdata)
        str(sg)
        numSigGenes(myGOdata)

        resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")


        allRes <- GenTable(myGOdata,
            Fisher = resultFisher,
            orderBy = "Fisher", ranksOf = "Fisher", topNodes = length(resultFisher@score)
        )
        allRes$Fisher <- as.numeric(allRes$Fisher)
        allRes$Term <- gsub(",", ";", allRes$Term)
        allRes$type <- GO_type_loop

        res_list[[GO_type_loop]] <- allRes[allRes$Fisher <= 0.01, ]
    }

    res_bind <- rbind(res_list[["BP"]], res_list[["CC"]], res_list[["MF"]])
    res_bind <- res_bind[!grepl("cellular_component|biological_process|molecular_function", res_bind$Term), ]

    return(res_bind)
}
### plot
gain_bind <- direction_fun(feture_res_table, "gain")
loss_bind <- direction_fun(feture_res_table, "loss")

################
### just BP
if (just_BP) {
    gain_bind = gain_bind %>% filter(type == "BP")
    loss_bind = loss_bind %>% filter(type == "BP")
}
file_name_addition = ifelse(just_BP, "_just_BP", "")
################

gain_col <- "#cf534c"
loss_col <- "#6397eb"

bubble_gain <- gain_bind %>%
    ggplot(aes(Significant, reorder(Term, Significant), size = Annotated, color = Fisher)) +
    scale_color_gradient("p.value", low = gain_col, high = "black") +
    labs(x = "Significant", y = "") +
    theme_bw() +
    theme(
        # plot.title=element_text(hjust=0.5),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 9.5),
        legend.position = "right",
        text = element_text(family = "serif")
    ) +
    geom_point() +
    facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
    guides(color = guide_colorbar(order = 1, barheight = 4))

bubble_loss <- loss_bind %>%
    ggplot(aes(Significant, reorder(Term, Significant), size = Annotated, color = Fisher)) +
    scale_color_gradient("p.value", low = loss_col, high = "black") + # theme_classic() +
    labs( # title = paste0(type_name," - ",gain_loss,"regulated transcripts"),
        x = "Significant", y = ""
    ) +
    theme_bw() +
    theme(
        # plot.title=element_text(hjust=0.5),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 9.5),
        legend.position = "right",
        text = element_text(family = "serif")
    ) +
    geom_point() +
    facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
    guides(color = guide_colorbar(order = 1, barheight = 4))

Height <- max(c(nrow(gain_bind), nrow(loss_bind))) / 6.25
if (Height < 3) {
    Height <- 3
}


source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/multiplot_ggplot2.R")

svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_correlation", file_name_addition, "_geneFeature.svg"), width = 9.90, height = Height, family = "serif")
multiplot(bubble_gain, bubble_loss, cols = 2)
dev.off()



##### revigo
# source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/revigo_function.R")
#
#      scatterPlot_fun <- Revigo_plots(
#          GO_df_up = res_bind,
#          GO_df_down = res_bind,
#          treatment = "mto1",
#          GO_type = "BP"
#      )
#
# svg(paste0(x.file, "up_ScatterPlot", x.file.suffix, ".svg"), width = 5, height = 2.38, family = "serif")
# print(x$scatterPlot_up)
# dev.off()
