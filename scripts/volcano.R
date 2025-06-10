plot_volcano <- function(res_obj,
                         FDR = 0.05,
                         xlim_both_sides = NULL,
                         ylim_up = NULL,
                         vlines = NULL,
                         title = NULL,
                         intgenes = NULL,
                         intgenes_color = "steelblue",
                         labels_intgenes = TRUE,
                         labels_repel = TRUE,
                         dot_size = 0.275,
                         alpha.yo = 0.4
                         ) {
  mydf <- as.data.frame(res_obj) %>% filter(!is.na(padj))
  mydf$id <- rownames(mydf)
  #mydf$isDE <- ifelse(is.na(res_obj$padj), FALSE, res_obj$padj < FDR)
  
  # Factor 'geneCat' with levels in the desired order
  mydf$geneCat <- with(mydf, ifelse(padj < FDR & log2FoldChange > 1, "Upregulated",
                                    ifelse(padj < FDR & log2FoldChange < -1, "Downregulated", "nonDE")))
  mydf$geneCat <- factor(mydf$geneCat, levels = c("Upregulated", "Downregulated", "nonDE"))
  
  mydf <- mydf[mydf$baseMean > 0, ]
  
  p <- ggplot(mydf, aes_string(x = "log2FoldChange", y = "-log10(padj)", color = "geneCat")) +
    geom_point(alpha = alpha.yo, size = dot_size)
  
  #p <- ggplot(mydf, aes_string(x = "log2FoldChange", y = "-log10(pvalue)")) +
  #  geom_point(aes_string(color = "isDE"), alpha = alpha.yo, size = 0.75)
  
  
  if (!is.null(ylim_up)) {
    p <- p + coord_cartesian(ylim = c(0, ylim_up), xlim = c(-xlim_both_sides, xlim_both_sides))
  } else {
    p <- p + coord_cartesian(ylim = c(0, 20), xlim = c(-xlim_both_sides, xlim_both_sides))
  }
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  p <- p + theme_classic() + #theme_bw() +
    scale_colour_manual(
      name = "Gene Expression", # This changes the legend title
      values = c("nonDE" = "gray60", "Upregulated" = "#a84848", "Downregulated" = "#5d60ba")
    ) +
    guides(color = guide_legend(override.aes = list(size = 2.5)))
    #scale_colour_manual(
    #  name = paste0("FDR = ", FDR),
    #  values = c("black", "#a84848"),
    #  labels = c("nonDE", "DE")
    #)
  
  y_intercept_value <- -log10(FDR)
  p <- p + 
    # Add a vertical line starting from y=5 at x=1
    geom_segment(aes(x = 1, y = y_intercept_value, xend = 1, yend = Inf), 
                 col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed") +
    # Add a horizontal line ending at x=-1 from the left
    geom_segment(aes(x = -1, y = y_intercept_value, xend = -1, yend = Inf), 
                 col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed") +
    # Add a horizontal line starting from x=1 to the right
    geom_segment(aes(x = 1, y = y_intercept_value, xend = Inf, yend = y_intercept_value), 
                 col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed") +
    geom_segment(aes(x = -Inf, y = y_intercept_value, xend = -1, yend = y_intercept_value), 
                 col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed")
  
  #p <- p + geom_vline(aes(xintercept = 1), col = "lightblue", alpha = 0.6, size = 1) +
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
  }
  
  p
}