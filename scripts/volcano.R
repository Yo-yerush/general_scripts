plot_volcano <- function(
    res_obj,
    FDR = 0.05,
    xlim_both_sides = NULL,
    ylim_up = NULL,
    vlines = NULL,
    title = NULL,
    y_lab = "-log10(padj)", 
    x_lab = "log2(FoldChange)",
    intgenes = NULL,
    intgenes_color = "steelblue",
    labels_intgenes = TRUE,
    labels_repel = TRUE,
    dot_size = 0.275,
    gene_groups = NULL, # groups: list vectors with IDs
    group_colors = NULL, # if one group, van use 'random' to not get alwayse green
    alpha.yo = ifelse(is.null(gene_groups), 0.4, 0.8)) {

  library(RColorBrewer)
  mydf <- as.data.frame(res_obj) %>% filter(!is.na(padj))
  mydf$id <- rownames(mydf)
  # mydf$isDE <- ifelse(is.na(res_obj$padj), FALSE, res_obj$padj < FDR)

  # Factor 'geneCat' with levels in the desired order
  ################## old
  #  mydf$geneCat <- with(mydf, ifelse(padj < FDR & log2FoldChange > 1, "Upregulated",
  #                                    ifelse(padj < FDR & log2FoldChange < -1, "Downregulated", "nonDE")))
  #  mydf$geneCat <- factor(mydf$geneCat, levels = c("Upregulated", "Downregulated", "nonDE"))
  ##################

  ################## new
  mydf$geneCat <- "nonDE"
  if (!is.null(gene_groups)) {
    for (grp in names(gene_groups)) {
      mydf$geneCat[mydf$gene_id %in% gene_groups[[grp]]] <- grp
    }
    lvl <- c(names(gene_groups), "nonDE")
    mydf$geneCat <- factor(mydf$geneCat, levels = lvl)
  } else {
    mydf$geneCat <- with(mydf, ifelse(padj < FDR & log2FoldChange > 1, "Upregulated",
      ifelse(padj < FDR & log2FoldChange < -1, "Downregulated", "nonDE")
    ))
    mydf$geneCat <- factor(mydf$geneCat, levels = c("Upregulated", "Downregulated", "nonDE"))
  }

  if (!is.null(gene_groups)) {
    if (length(gene_groups) == 1) {
      if (!is.null(group_colors)) {
        if (group_colors == "random") {
          color_pallet <- brewer.pal(8, "Set2")[-c(5,6,8)][sample(5, 1)]
        } else {
          color_pallet <- group_colors
        }
      } else {
         color_pallet <- brewer.pal(8, "Set2")
      }
    } else {
      color_set <- ifelse(length(gene_groups) > 4, "Set1", "Set2")
      if (length(gene_groups) < 5) {
        color_pallet <- brewer.pal(length(gene_groups), "Set2")
      } else if (length(gene_groups) < 9) {
        color_pallet <- brewer.pal(length(gene_groups), "Set1")
        if (!is.na(color_pallet[6])) {
          color_pallet[6] <- "#3b3b34"
        }
      } else {
        color_pallet <- rainbow(length(gene_groups))
      }
    }
    base_cols <- setNames(color_pallet, names(gene_groups))

    palette_vals <- c(base_cols, "nonDE" = "#b9b9b9")
    legend_name <- "Gene Groups"

  } else {
    palette_vals <- c("nonDE" = "gray60", "Upregulated" = "#a84848", "Downregulated" = "#5d60ba")
    legend_name <- "Gene Expression"
  }
  ##################

  if ("baseMean" %in% names(mydf)) {
    mydf <- mydf[mydf$baseMean > 0, ]
  }

  ##################

  mydf$geneCat <- gsub("_", " ", mydf$geneCat)
  names(palette_vals) <- gsub("_", " ", names(palette_vals))

  ##################

  p <- ggplot(mydf, aes_string(x = "log2FoldChange", y = "-log10(padj)", color = "geneCat")) +
    geom_point(alpha = alpha.yo, size = dot_size)

  # p <- ggplot(mydf, aes_string(x = "log2FoldChange", y = "-log10(pvalue)")) +
  #  geom_point(aes_string(color = "isDE"), alpha = alpha.yo, size = 0.75)

  max_y <- max(-log10(mydf$padj))
  ylim_up <- ifelse(max_y < 20, max_y, 20)
  if (!is.null(xlim_both_sides)) {
    p <- p + coord_cartesian(ylim = c(0, ylim_up), xlim = c(-xlim_both_sides, xlim_both_sides))
  } else {
    p <- p + coord_cartesian(ylim = c(0, ylim_up))
  }

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  p <- p + theme_classic() + # theme_bw() +
    scale_colour_manual(name = legend_name, values = palette_vals, breaks = names(palette_vals)) +
    guides(color = guide_legend(override.aes = list(size = 2.5)))
  # scale_colour_manual(
  #  name = paste0("FDR = ", FDR),
  #  values = c("black", "#a84848"),
  #  labels = c("nonDE", "DE")
  # )

  y_intercept_value <- -log10(FDR)
  p <- p +
    # Add a vertical line starting from y=5 at x=1
    geom_segment(aes(x = 1, y = y_intercept_value, xend = 1, yend = Inf),
      col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed"
    ) +
    # Add a horizontal line ending at x=-1 from the left
    geom_segment(aes(x = -1, y = y_intercept_value, xend = -1, yend = Inf),
      col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed"
    ) +
    # Add a horizontal line starting from x=1 to the right
    geom_segment(aes(x = 1, y = y_intercept_value, xend = Inf, yend = y_intercept_value),
      col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed"
    ) +
    geom_segment(aes(x = -Inf, y = y_intercept_value, xend = -1, yend = y_intercept_value),
      col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed"
    )

  # p <- p + geom_vline(aes(xintercept = 1), col = "lightblue", alpha = 0.6, size = 1) +
  #  geom_vline(aes(xintercept = -1), col = "lightblue", alpha = 0.6, size = 1) +
  #  geom_hline(aes(yintercept = -log10(0.01)), col = "lightblue", alpha = 0.6, size = 1)

  if (!is.null(intgenes)) {
    if ("symbol" %in% colnames(mydf)) {
      # use the gene names
      df_intgenes <- mydf[mydf$symbol %in% intgenes, ]
      df_intgenes$myids <- df_intgenes$symbol
    } else {
      # use whatever is there as id
      df_intgenes <- mydf[rownames(mydf) %in% intgenes, ]
      df_intgenes$myids <- rownames(df_intgenes)
    }

    # df_intgenes <- mydf[mydf$symbol %in% intgenes,]

    p <- p + geom_point(data = df_intgenes, aes_string("log2(FoldChange)", "-log10(padj)"), color = intgenes_color, size = 4)

    if (labels_intgenes) {
      if (labels_repel) {
        p <- p + geom_text_repel(
          data = df_intgenes, aes_string("log2(FoldChange)", "-log10(padj)", label = "myids"),
          color = intgenes_color, size = 5
        )
      } else {
        p <- p + geom_text(
          data = df_intgenes, aes_string("log2(FoldChange)", "-log10(padj)", label = "myids"),
          color = intgenes_color, size = 5, hjust = 0.25, vjust = -0.75
        )
      }
    }

    p <- p + labs(x = x_lab_name, y = y_lab_name)
  }

  p
}
