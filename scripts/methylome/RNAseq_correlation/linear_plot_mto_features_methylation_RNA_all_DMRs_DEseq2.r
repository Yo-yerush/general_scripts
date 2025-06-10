###############
## DMRs gene list

library(plyr)
library(dplyr)
library(corrplot)
library(ggpubr)
library(ggplot2)
library(lmerTest)
library(tidyr)
source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/multiplot_ggplot2.R")
source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/yo_theme_base_ggplot2.R")

gene_list_name = "stress_related_genes"
gene_list <- data.frame(gene_id = c(NULL))

# open new plot for some reason..........
ggplot() + theme_void()

is.transcript = F

annotation_type_names = c("promoter", "CDS", "intron", "five_prime_UTR", "three_prime_UTR", "transposable_element_gene")

#### addit annotation names ####
edit_ann_names <- function(x) {
    x <- gsub("^promoter$", "Promoters", x)
    x <- gsub("^gene$", "Genes", x)
    x <- gsub("^intron$", "Introns", x)
    x <- gsub("^five_prime_UTR$", "5'UTRs", x)
    x <- gsub("^three_prime_UTR$", "3'UTRs", x)
    x <- gsub("^transposable_element_gene$", "TEGs", x)
    return(x)
}
edit_ann_names_2 <- function(x) {
    x <- gsub("^promoter$", "Promoters", x)
    x <- gsub("^gene$", "Genes", x)
    x <- gsub("^intron$", "Introns", x)
    x <- gsub("^five_prime_UTR$", "fiveUTRs", x)
    x <- gsub("^three_prime_UTR$", "threeUTRs", x)
    x <- gsub("^transposable_element_gene$", "TEG", x)
    return(x)
}

# set list for all plots
annotation_plots_list <- setNames(vector("list", length(annotation_type_names)), annotation_type_names)
residuals_plots_list <- setNames(vector("list", length(annotation_type_names)), annotation_type_names)


for (treatment in c("mto1","mto3")) {
  for (context in c("CG","CHG","CHH")) {
    for (annotation_type in annotation_type_names) {
      path_2_save.0 = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/DMRs_genes/Linear_correlation/"
      path_2_save.1 = paste0(path_2_save.0, treatment, "/")
      path_2_save.2 = paste0(path_2_save.1, context, "/")

            # gene_list directory
            if (!is.null(gene_list$gene_id)) {
              path_2_save.2 = paste0(path_2_save.2, "/", gene_list_name, "/")
              dir.create(, showWarnings = F)
            }

      dir.create(path_2_save.0, showWarnings = F)
      dir.create(path_2_save.1, showWarnings = F)
      dir.create(path_2_save.2, showWarnings = F)

      if (treatment == "mto1") {
        mto_rnaseq_names = c("met14", "met15", "met16")
      } else {
        mto_rnaseq_names = c("met17", "met18", "met19")
      }
      wt_rnaseq_names = c("met20", "met22")

      # Load the methylation files of the promoters of interestL
      meth_matrix = read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/average.meth.DMR.levels/", treatment, "_vs_wt/", context, "/meth.", edit_ann_names_2(annotation_type), ".", context, ".", treatment, "_vs_wt.csv"))

      # Load the RNAseq '*.gene.results' files (output from RSEM pipeline)
      RNA.0 = read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/", treatment, "_vs_wt/norm.", treatment, "_vs_wt.DE.csv"))

      if (is.transcript) {
        names(RNA.0)[1] = "transcript_id"
        refseq_2_tair = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/TAIR10.1/TAIR_2_refseq_transcript.csv")[, 1:2]
        RNA = merge.data.frame(refseq_2_tair, RNA.0, by = "transcript_id")
      } else {
        names(RNA.0)[1] = "gene_id"
        RNA = RNA.0
      }


      #### Gene expression ####
      # edit table and names
      mto_names_pos = grep(paste(mto_rnaseq_names, collapse = "|"), names(RNA))
      names(RNA)[mto_names_pos] = paste0(treatment, ".", 1:3, "_RNA")

      wt_names_pos = grep(paste(wt_rnaseq_names, collapse = "|"), names(RNA))
      names(RNA)[wt_names_pos] = paste0("wt.", 1:2, "_RNA")

      RNA = RNA[, c(1, wt_names_pos, mto_names_pos)]


      #### Methylation ####
      meth_matrix = meth_matrix[, grep("gene_id|mto|wt", names(meth_matrix))]
      names(meth_matrix)[-1] = paste0(names(meth_matrix)[-1], "_meth")
      meth_matrix = na.omit(meth_matrix) # remove rows contain 'NA'

      #### Merge methylation and RNA dataframes ####
      meth_matrix_RNA = merge(meth_matrix, RNA, by = "gene_id", all = FALSE)

      # filter by gene_list (if exist)
      if (!is.null(gene_list$gene_id)) {
         meth_matrix_RNA = filter(meth_matrix_RNA, gene_id %in% gene_list$gene_id)
      }
      # meth_matrix_RNA <- meth_matrix_RNA %>%
      #  dplyr::select(gene_id, transcript_id, everything())

      meth_new = meth_matrix_RNA[, c(grep("gene_id", names(meth_matrix_RNA)), grep("meth", names(meth_matrix_RNA)))]
      RNA_new = meth_matrix_RNA[, c(grep("gene_id", names(meth_matrix_RNA)), grep("RNA", names(meth_matrix_RNA)))]

      # if (all(meth_new$gene_id == RNA_new$gene_id)) {
      #  merge_new = cbind(RNA_new[,-1:2], meth_new[,-(1:2)]) # for cor_df to save
      # }

      names(meth_new) = gsub("_meth", "", names(meth_new))
      names(RNA_new) = gsub("_RNA", "", names(RNA_new))

      # move the columns/rows:
      fin_meth = meth_new %>%
        pivot_longer(!c(gene_id), names_to = "Sample", values_to = "Meth")

      fin_RNA = RNA_new %>%
        pivot_longer(!c(gene_id), names_to = "Sample", values_to = "TPM")

      # check if both data frames are similar (for cbind)
      # if (
      #  all(fin_meth$Sample == fin_RNA$Sample) & all(fin_meth$gene_id == fin_RNA$gene_id)
      # ) {

      fin_meth_RNA = cbind(fin_meth, TPM = fin_RNA$TPM)
      # fin_meth_RNA <- fin_meth_RNA %>%
      #  dplyr::select(gene_id, everything())
      fin_meth_RNA$Meth[is.na(fin_meth_RNA$Meth)] = 0 # i dont know if its the right way

      #### final DF for plotting ####
      plot_df <- fin_meth_RNA %>%
        mutate(genotype = ifelse(grepl("^wt", Sample), "wt", treatment))

      # Add 1 constant to TPM values to avoid log(0)
      plot_df$TPM <- log(plot_df$TPM + 1)
      # plot_df$Meth <- log(plot_df$Meth + 1)

      #### linear model ####
      lm_model_0 = glm(TPM ~ Meth + genotype, data = plot_df)
      lm_model <- summary(lm_model_0)
      
      #### linear model p-Value ####
      pval_fun <- function(x) {
        ifelse(x <= 0.001, "***",
          ifelse(x <= 0.01, "**",
            ifelse(x <= 0.05, "*",
              "nf"
            )
          )
        )
      }
      M_pval = pval_fun(lm_model$coefficients["Meth", 4])
      G_pval = pval_fun(lm_model$coefficients["genotypewt", 4])

      #### linear model - fitted vs. residuals ####
      plot_df$residuals <- residuals(glm(TPM ~ Meth + genotype, data = plot_df))
      plot_df$fitted <- fitted(glm(TPM ~ Meth + genotype, data = plot_df))

      #### plot ####
      annotation_plots_list[[annotation_type]] = ggplot(plot_df, aes(x = Meth, y = TPM, color = genotype)) +
        geom_point(alpha = 0.25, size = 0.25) +
        #ggrastr::geom_point_rast(alpha = 0.25, size = 0.25) +
        stat_smooth(method = "glm", formula = y ~ x) +
        #stat_smooth() +
        scale_color_manual(values = c("#bf6828", "gray50")) +
        # 'M' and 'G' text
        annotate("text", x = Inf, y = Inf, label = "M", hjust = 4, vjust = 2, size = 3.5) +
        annotate("text", x = Inf, y = Inf, label = "G", hjust = 4.75, vjust = 4, size = 3.5) +
        # pValue text
        annotate("text", x = Inf, y = Inf, label = paste0("'", M_pval, "'"), hjust = 1.15, vjust = 2, size = 3.5) +
        annotate("text", x = Inf, y = Inf, label = paste0("'", G_pval, "'"), hjust = 1.15, vjust = 4, size = 3.5) +
        yo_theme_base() +
        theme(legend.position = "none") +
        labs(
          title = edit_ann_names(annotation_type),
          x = "Methylation levels",
          y = "log(Norm.Counts + 1)"
        )
      
      #### residuals plot ####
      residuals_plots_list[[annotation_type]] = ggplot(plot_df, aes(x = fitted, y = residuals, color = genotype)) +
        #ggrastr::geom_point_rast(alpha = 0.25, size = 0.25) +
        geom_point(alpha = 0.25, size = 0.25) +
        geom_smooth() +
        scale_color_manual(values = c("#bf6828", "gray50")) +
        yo_theme_base() +
        theme(legend.position = "none") +
        labs(
          title = edit_ann_names(annotation_type),
          x = "Fitted",
          y = "Residuals"
        )
      # write(lm_model, paste0(path_2_save.2, edit_ann_names(annotation_type), ".lm.stats.", context, ".", treatment, ".txt"))
      # }
    }
    
    
    #### legend ####
    p <- ggplot(plot_df, aes(x = Meth, y = TPM, fill = genotype)) +
      geom_point(alpha = 1, size = 6.5, shape = 21, stroke = 0.5, color = "black") +
      scale_fill_manual(
        values = c("mto1" = "#bf6828", "wt" = "gray50"),
        breaks = c(treatment, "wt") #, labels = c(expression(italic("mto1")), "wt")
        ) +
      ggthemes::theme_base() +
      theme(
        legend.text = element_text(size = 16),
        legend.position = "left"
      ) + # Increase legend text size and position it at the top left
      labs(fill = "")
    
    legend_p_0 <- ggplotGrob(p)$grobs[[which(sapply(ggplotGrob(p)$grobs, function(x) x$name) == "guide-box")]]
    legend_p <- ggplot() +
      theme_void() +
      annotation_custom(legend_p_0, xmin = 0, xmax = 0.5, ymin = 0.5, ymax = 1)
    # legend_p = grid::grid.draw(legend_p)

    #### save all plot as one, for each context ####

    ### main plot
    svg(paste0(path_2_save.2, "lm.stats.plot.", context, ".", treatment, ".svg"), width = 12, height = 5, family = "serif")
    multiplot(
      annotation_plots_list[[1]],
      annotation_plots_list[[2]],
      annotation_plots_list[[3]],
      annotation_plots_list[[4]],
      annotation_plots_list[[5]],
      annotation_plots_list[[6]],
      #annotation_plots_list[[7]],
      legend_p,
      cols = 4
    )
    dev.off()

    ### residuals plot
    svg(paste0(path_2_save.2, "lm.residuals.plot.", context, ".", treatment, ".svg"), width = 12, height = 5, family = "serif")
    multiplot(
      residuals_plots_list[[1]],
      residuals_plots_list[[2]],
      residuals_plots_list[[3]],
      residuals_plots_list[[4]],
      residuals_plots_list[[5]],
      residuals_plots_list[[6]],
      #residuals_plots_list[[7]],
      legend_p,
      cols = 4
    )
    dev.off()
  }
}