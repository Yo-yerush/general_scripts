library(dplyr)

# Define the treatments and contexts
treatments <- c("mto1", "mto3")
contexts <- c("CG", "CHG", "CHH")
annotations <- c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs", "TEG")

# Define colors for Positive and Negative correlations
correlation_colors <- c("Negative" = "#6c96d9", "Positive" = "#d96c6c")

# Set up the paths
base_path <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/"
output_path <- paste0(base_path, "NGS_merged_results/corr_with_methylations/corr_percentage_barplots/Gene_feature/")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Loop through each treatment
for (treatment in treatments) {
    # Create a subdirectory for each treatment
    treatment_path <- paste0(output_path, treatment, "/")
    dir.create(treatment_path, showWarnings = FALSE)

    ################
    # data frame for DEGs summerize table
    summerize_DEGs_with_cor = data.frame(
        annotation = annotations,
        positive_CG = NA, positive_CHG = NA, positive_CHH = NA,
        negative_CG = NA, negative_CHG = NA, negative_CHH = NA
    )
    ################

    # Loop through each context
    for (context in contexts) {
        # Initialize vectors to store percentage data for each annotation
        percentages_negative <- numeric(length(annotations))
        percentages_positive <- numeric(length(annotations))
        gene_counts_negative <- numeric(length(annotations))
        gene_counts_positive <- numeric(length(annotations))
        percentages_hyper_DMRs <- numeric(length(annotations))
        percentages_hyper_DMRs_negative <- numeric(length(annotations))
        percentages_hyper_DMRs_positive <- numeric(length(annotations))
        

        # Loop through each annotation to collect data
        for (i in seq_along(annotations)) {
            ann <- annotations[i]

            # Construct file paths for gene annotations
            annotation_file <- read.csv(paste0(
                base_path, "methylome_23/BSseq_results/", treatment, "_vs_wt/genome_annotation/",
                context, "/", ann, "_", context, "_genom_annotations.csv"
            ))

            # Read and filter gene IDs
            gene_ids <- annotation_file$gene_id
            filtered_gene_ids <- unique(gene_ids[grep("^AT[0-5]G", gene_ids)])
            total_genes <- length(filtered_gene_ids)

            # Construct file path for correlation data
            corr_file <- paste0(
                base_path, "NGS_merged_results/corr_with_methylations/by_DEseq2/Gene_feature/",
                treatment, "/", context, "/", ann, ".corr.", context, ".", treatment, ".csv"
            )

            # Read and filter correlation data
            corr_data <- read.csv(corr_file)
            corr_filtered <- subset(corr_data, pval < 0.05)

            # Calculate unique gene counts for negative and positive correlations
            neg_genes_list = data.frame(gene_id = unique(corr_filtered$gene_id[corr_filtered$cor < 0]))
            pos_genes_list = data.frame(gene_id = unique(corr_filtered$gene_id[corr_filtered$cor > 0]))
            neg_genes <- nrow(neg_genes_list)
            pos_genes <- nrow(pos_genes_list)
            
            # Calculate percentages from total
            percentages_negative[i] <- (neg_genes / total_genes) * 100
            percentages_positive[i] <- (pos_genes / total_genes) * 100

            # Store gene counts for annotation
            gene_counts_negative[i] <- neg_genes
            gene_counts_positive[i] <- pos_genes

            # hyper-DMRs percentage in positive/negative list
            ann_cor_merged = merge.data.frame(annotation_file, corr_filtered, by = "gene_id")
            pos_cor_df <- ann_cor_merged[ann_cor_merged$cor > 0, ]
            neg_cor_df = ann_cor_merged[ann_cor_merged$cor < 0, ]
            pos_cor_gain_count <- nrow(pos_cor_df[pos_cor_df$regionType == "gain", ])
            neg_cor_gain_count <- nrow(neg_cor_df[neg_cor_df$regionType == "gain", ])
            percentages_hyper_DMRs_positive[i] <- (pos_cor_gain_count / nrow(pos_cor_df)) * 100
            percentages_hyper_DMRs_negative[i] <- (neg_cor_gain_count / nrow(neg_cor_df)) * 100

            # hyper-DMRs percentage
            total_count <- nrow(ann_cor_merged)
            gain_count <- nrow(ann_cor_merged[ann_cor_merged$regionType == "gain", ])
            percentages_hyper_DMRs[i] <- (gain_count / total_count) * 100
            
            #########################
            # DEGs file (for DEGs 'summerize table')
            DEGs <- read.csv(paste0(
                base_path, "rnaseq_23/met23/", treatment, "_vs_wt/all.transcripts.", treatment, "_vs_wt.DE.csv"
            ))[, 1:4] %>%
                dplyr::rename(gene_id = "locus_tag") %>%
                dplyr::filter(pValue < 0.05) %>%
                dplyr::select(gene_id, log2FoldChange)
            
            neg_cor_n_DEGs = merge(DEGs, neg_genes_list, by = "gene_id")
            pos_cor_n_DEGs = merge(DEGs, pos_genes_list, by = "gene_id")

            neg_DEGs_total = neg_cor_n_DEGs %>% nrow()
            pos_DEGs_total = pos_cor_n_DEGs %>% nrow()

            neg_DEGs_downregulated = neg_cor_n_DEGs[neg_cor_n_DEGs$log2FoldChange < 0, ] %>% nrow()
            pos_DEGs_downregulated = pos_cor_n_DEGs[pos_cor_n_DEGs$log2FoldChange < 0, ] %>% nrow()

            neg_column = grep(paste0("negative_", context), names(summerize_DEGs_with_cor))
            pos_column = grep(paste0("positive_", context), names(summerize_DEGs_with_cor))
            
            neg_percentage = round((neg_DEGs_downregulated / neg_DEGs_total) * 100, 1)
            pos_percentage = round((pos_DEGs_downregulated / pos_DEGs_total) * 100, 1)

            summerize_DEGs_with_cor[i, neg_column] = paste0(neg_DEGs_downregulated, "/", neg_DEGs_total, " (", neg_percentage, "%)")
            summerize_DEGs_with_cor[i, pos_column] = paste0(pos_DEGs_downregulated, "/", pos_DEGs_total, " (", pos_percentage, "%)")

            summerize_DEGs_with_cor[i, neg_column] <- gsub("NaN%", "0%", summerize_DEGs_with_cor[i, neg_column])
            summerize_DEGs_with_cor[i, pos_column] <- gsub("NaN%", "0%", summerize_DEGs_with_cor[i, pos_column])
            #########################
        }

        # Create a data frame for plotting
        plot_data <- data.frame(
            Annotation = annotations,
            Negative = percentages_negative,
            Positive = percentages_positive
        )
        plot_data$Annotation = gsub("five", "5'", plot_data$Annotation)
        plot_data$Annotation = gsub("three", "3'", plot_data$Annotation)
        plot_data$Annotation = gsub("TEG", "TEGs", plot_data$Annotation)

        # 2-row "table" above the bars
        table_data <- data.frame(
            Row_pos = round(percentages_hyper_DMRs_positive, 0),
            Row_neg = round(percentages_hyper_DMRs_negative, 0)
        )
        table_data$Row_pos = gsub(NaN, "-", table_data$Row_pos)
        table_data$Row_neg = gsub(NaN, "-", table_data$Row_neg)

        # Determine the y-axis limit (max percentage + buffer)
        y_max <- max(plot_data$Negative + plot_data$Positive) * 1.1
        y_max <- ifelse(y_max > 100, 100, y_max) # Cap at 100%

        # Define bar positions on x-axis
        bar_positions = seq(1, 8.5, length.out = 6)
        dot_positions = seq(1, 9, length.out = 6)

        # Define the output SVG file name
        svg_filename <- paste0(treatment_path, treatment, "_", context, "_corr_percentage.svg")

        # Open SVG device
        svg(svg_filename, width = 3.25, height = 3.75, family = "serif")
        par(mar = c(4, 4, 3, 3.25))

        ### Create the barplot
        par(mgp = c(3, 0.5, 0))  
        barplot(
            height = t(as.matrix(plot_data[, c("Negative", "Positive")])),
            beside = FALSE,
            col = correlation_colors,
            border = "black",
            axes = FALSE,
            space = 0.5,
            ylim = c(0, 20),
            ylab = "",
            # main = paste0(context, " Context"),
            names.arg = plot_data$Annotation,
            cex.names = 0.8,
            las = 2
        )
        par(mgp = c(3, 1, 0))
        axis(side = 2, at = seq(0, 20, by = 5))
        mtext("Correlations (%)", side = 2, line = 3)

        ### correlated gene count in bars
        for (i in seq_along(annotations)) {
            text(
                x = bar_positions[i],
                y = plot_data$Negative[i],
                labels = ifelse(gene_counts_negative[i] == 0, "", gene_counts_negative[i]),
                pos = 1,
                cex = 0.7
            )
            is_pos_na <- is.na(plot_data$Positive[i] / (plot_data$Negative[i] + plot_data$Positive[i]))
            new_pos <- if (!is_pos_na && plot_data$Positive[i] / (plot_data$Negative[i] + plot_data$Positive[i]) > 0.1) 1 else 3
            text(
                x = bar_positions[i],
                y = plot_data$Negative[i] + plot_data$Positive[i],
                labels = ifelse(gene_counts_positive[i] == 0, "", gene_counts_positive[i]),
                pos = new_pos,
                cex = 0.7
            )
        }


        ### 2-row "table" above the bars
        par(xpd = TRUE) # Turn off clipping to allow drawing outside the plot
        table_base_y <- 21 # y-offset above the top of the y-axis (20)
        row_step <- 1.5 # distance between the two text rows

        for (i in seq_len(nrow(plot_data))) {
            # Row 1
            text(
                x = bar_positions[i],
                y = table_base_y,
                labels = table_data$Row_pos[i],
                cex = 0.7,
                pos = 3,
                col = correlation_colors[2],
                font = 2,
                offset = 0.2
            )

            # Row 2 (slightly above or below row 1)
            text(
                x = bar_positions[i],
                # y = table_base_y - row_step,
                y = table_base_y,
                labels = table_data$Row_neg[i],
                cex = 0.7,
                pos = 1,
                col = correlation_colors[1],
                font = 2,
                offset = 0.2
            )
        }

        # Draw a line between row 1 and row 2
        segments(
            x0 = min(bar_positions) - 0.5,
            y0 = table_base_y,
            x1 = max(bar_positions) + 0.5,
            y1 = table_base_y,
            col = "black",
            lwd = 0.4
        )

        # Draw a line between columns
        segments(
            x0 = bar_positions[1:5] + 0.75,
            y0 = table_base_y - 1,
            x1 = bar_positions[1:5] + 0.75,
            y1 = table_base_y + 1,
            col = "black",
            lwd = 0.4
        )

        # main title
        text(
            x = 0,
            y = table_base_y + 3,
            labels = paste0(context," context"),
            cex = 1,
            pos = 4,
            col = "black",
            font = 2
        )

        # hyper-DMRs legend for table
        text(
            x = max(bar_positions) + 0.4,
            y = table_base_y * 0.99,
            labels = "Hyper-\nDMRs (%)",
            cex = 0.75,
            pos = 4,
            col = "black",
            font = 2
        )

        # correlation legend
        text(
            x = max(bar_positions) + 0.4,
            y = table_base_y * 0.8,
            labels = "Correlation",
            cex = 0.75,
            pos = 4,
            col = "black",
            font = 2
        )
        
        legend(
            x = max(bar_positions) + 0.5,
            y = table_base_y * 0.85,
            legend = c("Positive", "Negative"),
            fill = rev(correlation_colors),
            border = "black",
            cex = 0.75,
            bty = "n",
            #title = expression(bold(" Correlation")),
            title = "",
            x.intersp = 0.5,
            xjust = 0,
        )

        par(xpd = FALSE)


        ### second y-axis data (hyper-DMRs in percentages)
        if (F) {
            par(new = TRUE)
            par(mar = c(5, 4, 3, 4))
            plot(
                x = dot_positions,
                y = percentages_hyper_DMRs,
                space = 0.5,
                type = "b",
                pch = 16,
                col = "#1f1e1ea2",
                axes = FALSE,
                xlab = "",
                ylab = "",
                ylim = c(0, 100),
                xlim = c(0.5, 9.5)
            )

            # barplot for the second y-axis values
            par(new = TRUE)
            barplot(
                numeric(length(annotations)),
                ylim = c(0, 100),
                ylab = "",
                col = NA,
                border = NA,
                axes = FALSE,
                names.arg = rep("", length(annotations))
            )
            axis(side = 4, at = seq(0, 100, by = 20))
            mtext("Hyper DMRs (%)", side = 4, line = 3)
        }

        # Close the SVG device
        dev.off()
    }

    #########################################
    ### DEGs summerize table
    write.csv(summerize_DEGs_with_cor, paste0(treatment_path, "cor_DEGs_summerize_table.csv"), row.names = F)
    #########################################
}

