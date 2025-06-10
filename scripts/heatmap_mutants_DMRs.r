library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(tibble)

options(width = 80)

dir.create("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/other_mutants_results/test_outputs", showWarnings = FALSE)
setwd("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/other_mutants_results/test_outputs")

# context <- "CHG"
# dmr_annotation <- "Genes"
# dmr_annotation <- "Transposable_Elements"

for (dmr_annotation in c("Genes", "Promoters", "Transposable_Elements")) {
    for (context in c("CG", "CHG", "CHH")) {
        load_mto <- function(x) {
            x_vs <- ifelse(grepl("mto", x), "_vs_wt", "_vs_EV")
            read.csv(paste0("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/", x, x_vs, "/genome_annotation/", context, "/", dmr_annotation, "_", context, "_genom_annotations.csv")) %>%
                select(gene_id, log2FC) %>%
                group_by(gene_id) %>%
                summarise(!!rlang::sym(x) := mean(log2FC, na.rm = TRUE)) %>%
                as.data.frame()
        }

        load_mutant <- function(x) {
            read.csv(paste0("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/other_mutants_results/Methylome.At_res/", x, "_vs_wt/genome_annotation/", context, "/", dmr_annotation, "_", context, "_genom_annotations.csv")) %>%
                select(gene_id, log2FC) %>%
                group_by(gene_id) %>%
                summarise(!!rlang::sym(x) := mean(log2FC, na.rm = TRUE)) %>%
                as.data.frame()
        }

        ##############################################################

        mto1_merged <- load_mto("mto1")

        for (i_name in c("mto3", "dCGS", "SSE_high", "SSE_low")) {
            mto1_merged <- merge(mto1_merged, load_mto(i_name), by = "gene_id", all.x = TRUE)
        }
        for (i_name in c("adcp1", "ddm1", "ago4", "dml3", "hen1", "suvh8", "cmt3", "dcl4", "nrpb2", "nrpd1", "nrpe1", "suvh4", "cmt2", "ibm1","met1", "ktf1")) { #
            mto1_merged <- merge(mto1_merged, load_mutant(i_name), by = "gene_id", all.x = TRUE)
        }


        # Prepare data for heatmap

        if (dmr_annotation == "Transposable_Elements") {
            TE_file <- read.csv(paste0("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Transposable_Elements_", context, "_genom_annotations.csv")) %>%
                select(gene_id, Transposon_Super_Family) %>%
                filter(grepl("Gypsy|Copia|LINE", Transposon_Super_Family)) %>%
                mutate(Transposon_Super_Family = gsub("LTR/|/L1|\\?", "", Transposon_Super_Family)) %>%
                distinct(gene_id, .keep_all = TRUE)

            mto1_merged <- merge(mto1_merged, TE_file, by = "gene_id")

            mto1_merged <- arrange(mto1_merged, Transposon_Super_Family, gene_id)
            # mto1_merged <- arrange(mto1_merged, Transposon_Super_Family, desc(mto1))

            heatmap_data <- mto1_merged[, !(names(mto1_merged) %in% c("gene_id", "Transposon_Super_Family"))]
            rownames(heatmap_data) <- mto1_merged$gene_id

            row_annotation <- data.frame(TEs = mto1_merged$Transposon_Super_Family)
            rownames(row_annotation) <- rownames(heatmap_data)
        } else {
            # mto1_merged <- arrange(mto1_merged, desc(mto1))
            mto1_merged <- arrange(mto1_merged, gene_id)

            heatmap_data <- mto1_merged[, !(names(mto1_merged) %in% c("gene_id"))]
            rownames(heatmap_data) <- mto1_merged$gene_id

            row_annotation <- NA
        }

        # Order columns by correlation with mto1
        # Calculate correlations, handling NA values properly
        mto1_cors <- numeric(ncol(heatmap_data))
        names(mto1_cors) <- colnames(heatmap_data)

        for (i in 1:ncol(heatmap_data)) {
            valid_pairs <- !is.na(heatmap_data[, "mto1"]) & !is.na(heatmap_data[, i]) &
                is.finite(heatmap_data[, "mto1"]) & is.finite(heatmap_data[, i])

            if (sum(valid_pairs) > 2) { # Need at least 3 valid pairs for correlation
                mto1_cors[i] <- cor(heatmap_data[valid_pairs, "mto1"], heatmap_data[valid_pairs, i], method = "pearson")
            } else {
                mto1_cors[i] <- NA # Not enough valid pairs
            }
        }

        col_order <- order(abs(mto1_cors), decreasing = TRUE, na.last = TRUE)
        heatmap_data <- heatmap_data[, col_order]

        svg(paste0(dmr_annotation, "_", context, "_heatmap_mutants_DMRs.svg"), width = 5, height = 8, family = "serif")

        pheatmap(as.matrix(heatmap_data),
            color = colorRampPalette(c("blue", "white", "red"))(100),
            breaks = seq(-1, 1, length.out = 101),
            scale = "none",
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            annotation_row = row_annotation,
            show_rownames = FALSE,
            main = paste0(dmr_annotation, " in ", context, " context"),
            na_col = "white"
        )
        dev.off()


        ##################################################################
        ##################################################################
        ##################################################################
        ##################################################################
        ##################################################################
        ##################################################################
        ##################################################################
        ##################################################################
        ##################################################################
        ##################################################################
        ##################################################################
        ##################################################################
        ###################### sonnet 3.7

        # Create scatter plots with regression lines to compare mto1 vs other columns
        plot_correlation_analysis <- function(heatmap_data) {
            # Create output directory
            plot_dir <- paste0("mto1_correlation_plots_", dmr_annotation, "_", context)
            if (!dir.exists(plot_dir)) {
                dir.create(plot_dir)
            }

            # Get column names excluding mto1
            mto1_col <- heatmap_data[, "mto1"]
            other_cols <- setdiff(colnames(heatmap_data), "mto1")

            # Create individual scatter plots (existing code)
            for (col_name in other_cols) {
                current_col_data <- heatmap_data[, col_name]

                # Filter out NA and Inf values for plotting and correlation
                valid_indices <- !is.na(mto1_col) & !is.na(current_col_data) &
                    is.finite(mto1_col) & is.finite(current_col_data)

                if (sum(valid_indices) > 5) { # Only proceed if we have enough data points
                    plot_data_mto1 <- mto1_col[valid_indices]
                    plot_data_current <- current_col_data[valid_indices]

                    # Calculate correlation
                    cor_value <- cor(plot_data_mto1, plot_data_current, method = "pearson")

                    # Fit linear model for regression line
                    lm_fit <- lm(plot_data_current ~ plot_data_mto1)

                    # Save plot
                    plot_file_name <- paste0(plot_dir, "/mto1_vs_", col_name, "_correlation.png")
                    png(plot_file_name, width = 1000, height = 800, res = 120, family = "serif")

                    # Create empty plot with appropriate limits
                    x_range <- range(plot_data_mto1, na.rm = TRUE)
                    y_range <- range(plot_data_current, na.rm = TRUE)

                    plot(plot_data_mto1, plot_data_current,
                        xlab = "mto1 log2FC",
                        ylab = paste(col_name, "log2FC"),
                        main = paste0(
                            "mto1 vs ", col_name, " (", dmr_annotation, " - ", context, ")\nr = ", round(cor_value, 3),
                            ", n = ", sum(valid_indices)
                        ),
                        pch = 16,
                        col = rgb(0, 0, 1, 0.4),
                        xlim = x_range,
                        ylim = y_range
                    )

                    # Add diagonal line (y=x) for reference
                    abline(a = 0, b = 1, col = "darkgray", lty = 2, lwd = 1.5)

                    # Add regression line
                    abline(lm_fit, col = "red", lwd = 2)

                    # Add grid and legend
                    grid()
                    legend("topleft",
                        legend = c("y = x", paste0(
                            "Regression: y = ", round(lm_fit$coefficients[1], 3), " + ",
                            round(lm_fit$coefficients[2], 3), "x"
                        )),
                        col = c("darkgray", "red"),
                        lty = c(2, 1),
                        lwd = c(1.5, 2),
                        bg = "white"
                    )

                    dev.off()
                    cat("Created correlation plot:", plot_file_name, "\n")
                } else {
                    cat("Skipping", col_name, "- not enough valid data points for correlation analysis\n")
                }
            }

            # Create all correlation summary plots
            create_correlation_matrix(heatmap_data, plot_dir)
            create_correlation_barplot(heatmap_data, plot_dir)
            create_correlation_categories_boxplot(heatmap_data, plot_dir)
            create_ranked_correlation_plot(heatmap_data, plot_dir)
        }

        # Create a correlation matrix heatmap for all columns
        create_correlation_matrix <- function(heatmap_data, plot_dir) {
            # Remove rows with NA or Inf values
            clean_data <- heatmap_data[complete.cases(heatmap_data), ]

            if (nrow(clean_data) >= 5) { # Only proceed if we have enough complete rows
                # Calculate correlation matrix
                cor_matrix <- cor(clean_data, method = "pearson")

                # Save correlation matrix plot
                cor_plot_file <- paste0(plot_dir, "/correlation_matrix_", dmr_annotation, "_", context, ".png")
                png(cor_plot_file, width = 1000, height = 900, res = 120, family = "serif")

                # Use pheatmap for better visualization
                pheatmap(cor_matrix,
                    color = colorRampPalette(c("blue", "white", "red"))(100),
                    breaks = seq(-1, 1, length.out = 101),
                    display_numbers = TRUE,
                    number_format = "%.2f",
                    fontsize_number = 8,
                    fontsize = 10,
                    main = paste("DMR Profile Correlation Matrix -", dmr_annotation, context)
                )

                dev.off()
                cat("Created correlation matrix plot:", cor_plot_file, "\n")
            } else {
                cat("Skipping correlation matrix - not enough complete rows\n")
            }
        }

        # Function to create density plots for each comparison
        create_density_plots <- function(heatmap_data) {
            # Create output directory
            plot_dir <- paste0("mto1_density_plots_", dmr_annotation, "_", context)
            if (!dir.exists(plot_dir)) {
                dir.create(plot_dir)
            }

            # Calculate log2FC differences between mto1 and others
            mto1_col <- heatmap_data[, "mto1"]
            other_cols <- setdiff(colnames(heatmap_data), "mto1")

            # Combined plot for all differences
            plot_file_name <- paste0(plot_dir, "/all_differences_density_", dmr_annotation, "_", context, ".png")
            png(plot_file_name, width = 1200, height = 800, res = 120, family = "serif")

            # Setup the plot
            plot(0, 0,
                type = "n", xlim = c(-2, 2), ylim = c(0, 2),
                xlab = "log2FC Difference (other - mto1)",
                ylab = "Density",
                main = paste("Distribution of log2FC Differences -", dmr_annotation, context)
            )

            colors <- rainbow(length(other_cols), alpha = 0.7)
            legend_text <- character()

            for (i in seq_along(other_cols)) {
                col_name <- other_cols[i]
                current_col_data <- heatmap_data[, col_name]

                # Calculate differences where both values are valid
                valid_indices <- !is.na(mto1_col) & !is.na(current_col_data) &
                    is.finite(mto1_col) & is.finite(current_col_data)

                if (sum(valid_indices) > 10) { # Only proceed if we have enough data points
                    differences <- current_col_data[valid_indices] - mto1_col[valid_indices]

                    # Add density curve
                    dens <- density(differences, na.rm = TRUE, from = -2, to = 2)
                    lines(dens, col = colors[i], lwd = 2)

                    # Add this to legend
                    legend_text[i] <- paste0(col_name, " (n=", sum(valid_indices), ")")
                } else {
                    legend_text[i] <- paste0(col_name, " (insufficient data)")
                }
            }

            # Add vertical line at 0 (no difference)
            abline(v = 0, lty = 2, col = "darkgray", lwd = 1.5)

            # Add legend
            legend("topright", legend = legend_text, col = colors, lwd = 2, cex = 0.8, bg = "white")

            dev.off()
            cat("Created density plot of differences:", plot_file_name, "\n")

            # Create boxplot of differences
            create_boxplot(heatmap_data, plot_dir)
        }

        # Create boxplot of differences
        create_boxplot <- function(heatmap_data, plot_dir) {
            mto1_col <- heatmap_data[, "mto1"]
            other_cols <- setdiff(colnames(heatmap_data), "mto1")

            differences_list <- list()
            labels <- character()

            for (col_name in other_cols) {
                current_col_data <- heatmap_data[, col_name]

                valid_indices <- !is.na(mto1_col) & !is.na(current_col_data) &
                    is.finite(mto1_col) & is.finite(current_col_data)

                if (sum(valid_indices) > 5) {
                    differences <- current_col_data[valid_indices] - mto1_col[valid_indices]
                    differences_list[[col_name]] <- differences
                    labels <- c(labels, col_name)
                }
            }

            if (length(differences_list) > 0) {
                plot_file_name <- paste0(plot_dir, "/differences_boxplot_", dmr_annotation, "_", context, ".png")
                png(plot_file_name, width = 1200, height = 800, res = 120, family = "serif")

                par(mar = c(8, 4, 4, 2))
                boxplot(differences_list,
                    main = paste("Differences in log2FC compared to mto1 -", dmr_annotation, context),
                    ylab = "log2FC Difference (other - mto1)",
                    col = rainbow(length(differences_list), alpha = 0.7),
                    las = 2
                ) # Rotate x-axis labels

                # Add horizontal line at 0
                abline(h = 0, lty = 2, col = "darkgray", lwd = 1.5)

                dev.off()
                cat("Created boxplot of differences:", plot_file_name, "\n")
            } else {
                cat("Skipping boxplot - insufficient data\n")
            }
        }

        # Create bar plot showing correlation coefficients with mto1
        create_correlation_barplot <- function(heatmap_data, plot_dir) {
            mto1_col <- heatmap_data[, "mto1"]
            other_cols <- setdiff(colnames(heatmap_data), "mto1")

            correlations <- numeric()
            sample_sizes <- numeric()
            col_names <- character()

            # Calculate correlations for each column with mto1
            for (col_name in other_cols) {
                current_col_data <- heatmap_data[, col_name]

                # Filter out NA and Inf values
                valid_indices <- !is.na(mto1_col) & !is.na(current_col_data) &
                    is.finite(mto1_col) & is.finite(current_col_data)

                if (sum(valid_indices) > 5) {
                    cor_value <- cor(mto1_col[valid_indices], current_col_data[valid_indices], method = "pearson")

                    # Check if correlation is valid (not NA or Inf)
                    if (!is.na(cor_value) && is.finite(cor_value)) {
                        correlations <- c(correlations, cor_value)
                        sample_sizes <- c(sample_sizes, sum(valid_indices))
                        col_names <- c(col_names, col_name)
                    } else {
                        cat("Warning: Invalid correlation for", col_name, "\n")
                    }
                }
            }

            if (length(correlations) > 0) {
                # Create data frame and sort by absolute correlation
                cor_df <- data.frame(
                    mutant = col_names,
                    correlation = correlations,
                    abs_correlation = abs(correlations),
                    sample_size = sample_sizes
                )
                cor_df <- cor_df[order(-cor_df$abs_correlation), ]

                # Create bar plot
                plot_file_name <- paste0(plot_dir, "/mto1_correlation_barplot_", dmr_annotation, "_", context, ".png")
                png(plot_file_name, width = 1200, height = 800, res = 120, family = "serif")

                par(mar = c(10, 4, 4, 2))

                # Color bars based on positive/negative correlation
                bar_colors <- ifelse(cor_df$correlation > 0, "steelblue", "orange")

                barplot(cor_df$correlation,
                    names.arg = paste0(cor_df$mutant, "\n(n=", cor_df$sample_size, ")"),
                    main = paste("Correlation with mto1 -", dmr_annotation, context),
                    ylab = "Pearson Correlation Coefficient",
                    col = bar_colors,
                    las = 2,
                    ylim = c(-1, 1)
                )

                # Add horizontal lines for reference
                abline(h = 0, lty = 1, col = "black", lwd = 1)
                abline(h = c(-0.5, 0.5), lty = 2, col = "gray", lwd = 1)

                # Add text labels on bars
                text(
                    x = 1:length(cor_df$correlation),
                    y = cor_df$correlation + ifelse(cor_df$correlation > 0, 0.05, -0.05),
                    labels = round(cor_df$correlation, 2),
                    cex = 0.8
                )

                # Add legend
                legend("topright",
                    legend = c("Positive correlation", "Negative correlation"),
                    fill = c("steelblue", "orange"),
                    bg = "white"
                )

                dev.off()
                cat("Created correlation bar plot:", plot_file_name, "\n")

                # Print correlation summary
                cat("\nCorrelation Summary (sorted by absolute correlation):\n")
                print(cor_df[, c("mutant", "correlation", "sample_size")])

                return(cor_df)
            } else {
                cat("No valid correlations to plot\n")
                return(NULL)
            }
        }

        # Fix the correlation matrix function as well
        create_correlation_matrix <- function(heatmap_data, plot_dir) {
            # Remove rows with NA or Inf values more carefully
            # First, replace Inf with NA
            clean_data <- heatmap_data
            clean_data[!is.finite(as.matrix(clean_data))] <- NA

            # Remove rows with all NA
            clean_data <- clean_data[rowSums(!is.na(clean_data)) > 0, ]

            # Only keep columns that have enough data
            col_keep <- colSums(!is.na(clean_data)) > 5
            clean_data <- clean_data[, col_keep, drop = FALSE]

            if (nrow(clean_data) >= 5 && ncol(clean_data) >= 2) {
                # Calculate correlation matrix with pairwise complete observations
                cor_matrix <- cor(clean_data, method = "pearson", use = "pairwise.complete.obs")

                # Check for any remaining NA/Inf values in correlation matrix
                if (any(!is.finite(cor_matrix))) {
                    cat("Warning: Non-finite values in correlation matrix, cleaning...\n")
                    cor_matrix[!is.finite(cor_matrix)] <- 0
                }

                # Save correlation matrix plot
                cor_plot_file <- paste0(plot_dir, "/correlation_matrix_", dmr_annotation, "_", context, ".png")
                png(cor_plot_file, width = 1000, height = 900, res = 120, family = "serif")

                # Use pheatmap for better visualization
                tryCatch(
                    {
                        pheatmap(cor_matrix,
                            color = colorRampPalette(c("blue", "white", "red"))(100),
                            breaks = seq(-1, 1, length.out = 101),
                            display_numbers = TRUE,
                            number_format = "%.2f",
                            fontsize_number = 8,
                            fontsize = 10,
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            main = paste("DMR Profile Correlation Matrix -", dmr_annotation, context)
                        )
                    },
                    error = function(e) {
                        cat("Error in pheatmap:", e$message, "\n")
                        cat("Creating simple correlation plot instead...\n")

                        # Fallback to base R plotting
                        image(1:ncol(cor_matrix), 1:nrow(cor_matrix), t(cor_matrix),
                            col = colorRampPalette(c("blue", "white", "red"))(100),
                            xlab = "", ylab = "", main = paste("DMR Profile Correlation Matrix -", dmr_annotation, context),
                            axes = FALSE
                        )
                        axis(1, at = 1:ncol(cor_matrix), labels = colnames(cor_matrix), las = 2)
                        axis(2, at = 1:nrow(cor_matrix), labels = rownames(cor_matrix), las = 2)
                    }
                )

                dev.off()
                cat("Created correlation matrix plot:", cor_plot_file, "\n")
            } else {
                cat("Skipping correlation matrix - not enough complete data\n")
            }
        }

        # Create correlation comparison boxplot (grouped by correlation strength)
        create_correlation_categories_boxplot <- function(heatmap_data, plot_dir) {
            mto1_col <- heatmap_data[, "mto1"]
            other_cols <- setdiff(colnames(heatmap_data), "mto1")

            # Categorize mutants by correlation strength
            high_pos_cor <- list()
            moderate_pos_cor <- list()
            low_cor <- list()
            moderate_neg_cor <- list()
            high_neg_cor <- list()

            category_names <- c(
                "High Positive\n(r > 0.7)", "Moderate Positive\n(0.3 < r ≤ 0.7)",
                "Low\n(|r| ≤ 0.3)", "Moderate Negative\n(-0.7 ≤ r < -0.3)",
                "High Negative\n(r < -0.7)"
            )

            for (col_name in other_cols) {
                current_col_data <- heatmap_data[, col_name]

                valid_indices <- !is.na(mto1_col) & !is.na(current_col_data) &
                    is.finite(mto1_col) & is.finite(current_col_data)

                if (sum(valid_indices) > 5) {
                    cor_value <- cor(mto1_col[valid_indices], current_col_data[valid_indices], method = "pearson")

                    # Categorize based on correlation strength
                    if (cor_value > 0.7) {
                        high_pos_cor[[col_name]] <- cor_value
                    } else if (cor_value > 0.3) {
                        moderate_pos_cor[[col_name]] <- cor_value
                    } else if (abs(cor_value) <= 0.3) {
                        low_cor[[col_name]] <- cor_value
                    } else if (cor_value >= -0.7) {
                        moderate_neg_cor[[col_name]] <- cor_value
                    } else {
                        high_neg_cor[[col_name]] <- cor_value
                    }
                }
            }

            # Combine into list for boxplot
            correlation_categories <- list(
                unlist(high_pos_cor),
                unlist(moderate_pos_cor),
                unlist(low_cor),
                unlist(moderate_neg_cor),
                unlist(high_neg_cor)
            )

            # Remove empty categories
            non_empty <- sapply(correlation_categories, length) > 0
            correlation_categories <- correlation_categories[non_empty]
            category_names <- category_names[non_empty]

            if (length(correlation_categories) > 0) {
                plot_file_name <- paste0(plot_dir, "/correlation_categories_boxplot_", dmr_annotation, "_", context, ".png")
                png(plot_file_name, width = 1000, height = 800, res = 120, family = "serif")

                par(mar = c(8, 4, 4, 2))

                boxplot(correlation_categories,
                    names = category_names,
                    main = paste("Distribution of Correlations with mto1 by Category -", dmr_annotation, context),
                    ylab = "Correlation Coefficient",
                    col = c("darkgreen", "lightgreen", "gray", "orange", "red")[1:length(correlation_categories)],
                    las = 2
                )

                # Add horizontal line at 0
                abline(h = 0, lty = 2, col = "black", lwd = 1)

                # Add sample sizes to plot
                for (i in 1:length(correlation_categories)) {
                    text(i, par("usr")[3] - 0.1,
                        paste("n =", length(correlation_categories[[i]])),
                        cex = 0.8
                    )
                }

                dev.off()
                cat("Created correlation categories boxplot:", plot_file_name, "\n")
            } else {
                cat("No correlations found for category boxplot\n")
            }
        }

        # Create a ranked correlation plot (horizontal bar chart for better readability)
        create_ranked_correlation_plot <- function(heatmap_data, plot_dir) {
            mto1_col <- heatmap_data[, "mto1"]
            other_cols <- setdiff(colnames(heatmap_data), "mto1")

            correlations <- numeric()
            sample_sizes <- numeric()
            col_names <- character()

            for (col_name in other_cols) {
                current_col_data <- heatmap_data[, col_name]

                valid_indices <- !is.na(mto1_col) & !is.na(current_col_data) &
                    is.finite(mto1_col) & is.finite(current_col_data)

                if (sum(valid_indices) > 5) {
                    cor_value <- cor(mto1_col[valid_indices], current_col_data[valid_indices], method = "pearson")
                    correlations <- c(correlations, cor_value)
                    sample_sizes <- c(sample_sizes, sum(valid_indices))
                    col_names <- c(col_names, col_name)
                }
            }

            if (length(correlations) > 0) {
                # Sort by correlation value (not absolute)
                cor_df <- data.frame(
                    mutant = col_names,
                    correlation = correlations,
                    sample_size = sample_sizes
                )
                cor_df <- cor_df[order(cor_df$correlation), ]

                plot_file_name <- paste0(plot_dir, "/ranked_correlation_plot_", dmr_annotation, "_", context, ".png")
                png(plot_file_name, width = 1000, height = 800, res = 120, family = "serif")

                par(mar = c(4, 8, 4, 2))

                # Create horizontal bar plot
                bar_colors <- ifelse(cor_df$correlation > 0, "steelblue", "orange")

                barplot(cor_df$correlation,
                    names.arg = paste0(cor_df$mutant, " (n=", cor_df$sample_size, ")"),
                    main = paste("Ranked Correlations with mto1 -", dmr_annotation, context),
                    xlab = "Pearson Correlation Coefficient",
                    col = bar_colors,
                    horiz = TRUE,
                    las = 1,
                    xlim = c(-1, 1)
                )

                # Add vertical lines for reference
                abline(v = 0, lty = 1, col = "black", lwd = 2)
                abline(v = c(-0.5, 0.5), lty = 2, col = "gray", lwd = 1)
                abline(v = c(-0.3, 0.3), lty = 3, col = "lightgray", lwd = 1)

                dev.off()
                cat("Created ranked correlation plot:", plot_file_name, "\n")
            }
        }

        # Run all the plotting functions
        if (exists("heatmap_data") && ncol(heatmap_data) > 1) {
            plot_correlation_analysis(heatmap_data)
            create_density_plots(heatmap_data)

            cat("\nAll plots have been created successfully in the designated directories.\n")
        } else {
            cat("Error: heatmap_data not found or has insufficient columns for analysis.\n")
        }
    }
}

###############################################
###############################################
###############################################

# Create a better organized comprehensive summary plot
create_comprehensive_summary_plot <- function() {
  # Need to add required libraries
  library(tidyr)
  
  # Initialize data collection
  all_correlations <- data.frame()
  
  # Collect correlation data from all combinations
  for (dmr_annotation in c("Genes", "Promoters", "Transposable_Elements")) {
    for (context in c("CG", "CHG", "CHH")) {
      
      # Load and prepare data (same as in your main loop)
      load_mto <- function(x) {
        x_vs <- ifelse(grepl("mto", x), "_vs_wt", "_vs_EV")
        read.csv(paste0("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/", x, x_vs, "/genome_annotation/", context, "/", dmr_annotation, "_", context, "_genom_annotations.csv")) %>%
          select(gene_id, log2FC) %>%
          group_by(gene_id) %>%
          summarise(!!rlang::sym(x) := mean(log2FC, na.rm = TRUE)) %>%
          as.data.frame()
      }
      
      load_mutant <- function(x) {
        read.csv(paste0("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/other_mutants_results/Methylome.At_res/", x, "_vs_wt/genome_annotation/", context, "/", dmr_annotation, "_", context, "_genom_annotations.csv")) %>%
          select(gene_id, log2FC) %>%
          group_by(gene_id) %>%
          summarise(!!rlang::sym(x) := mean(log2FC, na.rm = TRUE)) %>%
          as.data.frame()
      }
      
      # Build merged dataset
      mto1_merged <- load_mto("mto1")
      
      for (i_name in c("mto3", "dCGS", "SSE_high", "SSE_low")) {
        mto1_merged <- merge(mto1_merged, load_mto(i_name), by = "gene_id", all.x = TRUE)
      }
      for (i_name in c("adcp1", "ddm1", "ago4", "dml3", "hen1", "suvh8", "cmt3", "dcl4", "nrpb2", "nrpd1", "nrpe1", "suvh4", "cmt2", "ibm1")) {
        mto1_merged <- merge(mto1_merged, load_mutant(i_name), by = "gene_id", all.x = TRUE)
      }
      
      # Prepare heatmap data
      heatmap_data <- mto1_merged[, !(names(mto1_merged) %in% c("gene_id"))]
      
      # Calculate correlations with mto1
      mto1_col <- heatmap_data[, "mto1"]
      other_cols <- setdiff(colnames(heatmap_data), "mto1")
      
      for (col_name in other_cols) {
        current_col_data <- heatmap_data[, col_name]
        
        valid_indices <- !is.na(mto1_col) & !is.na(current_col_data) & 
                        is.finite(mto1_col) & is.finite(current_col_data)
        
        if (sum(valid_indices) > 5) {
          cor_value <- cor(mto1_col[valid_indices], current_col_data[valid_indices], method = "pearson")
          
          if (!is.na(cor_value) && is.finite(cor_value)) {
            all_correlations <- rbind(all_correlations, data.frame(
              mutant = col_name,
              correlation = cor_value,
              context = context,
              annotation = dmr_annotation,
              sample_size = sum(valid_indices),
              condition = paste(dmr_annotation, context, sep = "_"),
              context_group = ifelse(context == "CG", "CG", "non-CG"),
              annotation_group = case_when(
                dmr_annotation == "Genes" ~ "Genes",
                dmr_annotation == "Promoters" ~ "Promoters", 
                dmr_annotation == "Transposable_Elements" ~ "TEs"
              )
            ))
          }
        }
      }
    }
  }
  
  # Create three separate comprehensive plots
  
  # 1. CG Context Summary Plot
  png("CG_CONTEXT_summary_mto1_correlations.png", 
      width = 1600, height = 1000, res = 120, family = "serif")
  
  layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE))
  
  # Filter CG data
  cg_data <- all_correlations[all_correlations$context == "CG", ]
  
  # Panel 1: CG correlations by annotation
  par(mar = c(8, 4, 4, 2))
  cg_avg <- cg_data %>%
    group_by(mutant, annotation_group) %>%
    summarise(mean_cor = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(abs(mean_cor)))
  
  # Create grouped barplot
  cg_matrix <- cg_avg %>%
    pivot_wider(names_from = annotation_group, values_from = mean_cor, values_fill = 0) %>%
    column_to_rownames("mutant") %>%
    as.matrix()
  
  barplot(t(cg_matrix), beside = TRUE, 
          main = "CG Context: Correlations by Annotation",
          ylab = "Correlation with mto1", 
          col = c("lightblue", "lightgreen", "lightcoral"),
          legend.text = c("Genes", "Promoters", "TEs"),
          las = 2, ylim = c(-1, 1))
  abline(h = 0, lty = 1, col = "black")
  
  # Panel 2: CG correlation distribution
  par(mar = c(6, 4, 4, 2))
  cg_by_annotation <- split(cg_data$correlation, cg_data$annotation_group)
  boxplot(cg_by_annotation,
          main = "CG Context: Correlation Distribution",
          ylab = "Correlation with mto1",
          col = c("lightblue", "lightcoral", "lightgreen"),
          ylim = c(-1, 1))
  abline(h = 0, lty = 2, col = "darkgray")
  
  # Panel 3: Top CG correlations
  par(mar = c(4, 8, 4, 2))
  top_cg <- cg_data %>%
    group_by(mutant) %>%
    summarise(mean_cor = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(abs(mean_cor))) %>%
    head(10) %>%
    arrange(mean_cor)
  
  barplot(top_cg$mean_cor,
          names.arg = top_cg$mutant,
          main = "Top 10 CG Correlations",
          xlab = "Correlation",
          col = ifelse(top_cg$mean_cor > 0, "steelblue", "orange"),
          horiz = TRUE, las = 1, xlim = c(-1, 1))
  abline(v = 0, lty = 1, col = "black", lwd = 2)
  
  # Panel 4: CG heatmap
  par(mar = c(8, 8, 4, 2))
  cg_heatmap_data <- cg_data %>%
    select(mutant, correlation, annotation_group) %>%
    pivot_wider(names_from = annotation_group, values_from = correlation, values_fill = NA) %>%
    column_to_rownames("mutant") %>%
    as.matrix()
  
  image(1:ncol(cg_heatmap_data), 1:nrow(cg_heatmap_data),
        t(cg_heatmap_data[nrow(cg_heatmap_data):1, ]),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        main = "CG Context Heatmap", axes = FALSE, zlim = c(-1, 1))
  axis(1, at = 1:ncol(cg_heatmap_data), labels = colnames(cg_heatmap_data), las = 2)
  axis(2, at = 1:nrow(cg_heatmap_data), labels = rev(rownames(cg_heatmap_data)), las = 2, cex.axis = 0.7)
  
  dev.off()
  
  # 2. Non-CG Context Summary Plot  
  png("NON_CG_CONTEXT_summary_mto1_correlations.png", 
      width = 1600, height = 1200, res = 120, family = "serif")
  
  layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE))
  
  # Filter non-CG data
  non_cg_data <- all_correlations[all_correlations$context != "CG", ]
  
  # Panel 1: CHG correlations by annotation
  par(mar = c(8, 4, 4, 2))
  chg_data <- non_cg_data[non_cg_data$context == "CHG", ]
  chg_by_annotation <- split(chg_data$correlation, chg_data$annotation_group)
  boxplot(chg_by_annotation,
          main = "CHG Context by Annotation",
          ylab = "Correlation with mto1",
          col = c("lightblue", "lightcoral", "lightgreen"),
          las = 2, ylim = c(-1, 1))
  abline(h = 0, lty = 2, col = "darkgray")
  
  # Panel 2: CHH correlations by annotation  
  par(mar = c(8, 4, 4, 2))
  chh_data <- non_cg_data[non_cg_data$context == "CHH", ]
  chh_by_annotation <- split(chh_data$correlation, chh_data$annotation_group)
  boxplot(chh_by_annotation,
          main = "CHH Context by Annotation", 
          ylab = "Correlation with mto1",
          col = c("lightblue", "lightcoral", "lightgreen"),
          las = 2, ylim = c(-1, 1))
  abline(h = 0, lty = 2, col = "darkgray")
  
  # Panel 3: Non-CG average correlations
  par(mar = c(8, 4, 4, 2))
  non_cg_avg <- non_cg_data %>%
    group_by(mutant) %>%
    summarise(mean_cor = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(abs(mean_cor))) %>%
    head(15)
  
  barplot(non_cg_avg$mean_cor,
          names.arg = non_cg_avg$mutant,
          main = "Non-CG Average Correlations",
          ylab = "Mean Correlation",
          col = ifelse(non_cg_avg$mean_cor > 0, "steelblue", "orange"),
          las = 2, ylim = c(-1, 1))
  abline(h = 0, lty = 1, col = "black")
  
  # Panel 4: CHG vs CHH comparison
  par(mar = c(6, 4, 4, 2))
  context_comparison <- split(non_cg_data$correlation, non_cg_data$context)
  boxplot(context_comparison,
          main = "CHG vs CHH Comparison",
          ylab = "Correlation with mto1", 
          col = c("lightpink", "lightcyan"),
          ylim = c(-1, 1))
  abline(h = 0, lty = 2, col = "darkgray")
  
  # Panel 5: Top non-CG by context
  par(mar = c(4, 8, 4, 2))
  top_non_cg_chg <- chg_data %>%
    group_by(mutant) %>%
    summarise(mean_cor = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(abs(mean_cor))) %>%
    head(8)
  
  barplot(top_non_cg_chg$mean_cor,
          names.arg = top_non_cg_chg$mutant,
          main = "Top CHG Correlations",
          xlab = "Correlation",
          col = ifelse(top_non_cg_chg$mean_cor > 0, "steelblue", "orange"),
          horiz = TRUE, las = 1, xlim = c(-1, 1))
  abline(v = 0, lty = 1, col = "black", lwd = 2)
  
  # Panel 6: Non-CG heatmap by context
  par(mar = c(8, 8, 4, 2))
  non_cg_heatmap <- non_cg_data %>%
    select(mutant, correlation, context, annotation_group) %>%
    unite("condition", context, annotation_group, sep = "_") %>%
    pivot_wider(names_from = condition, values_from = correlation, values_fill = NA) %>%
    column_to_rownames("mutant") %>%
    as.matrix()
  
  image(1:ncol(non_cg_heatmap), 1:nrow(non_cg_heatmap),
        t(non_cg_heatmap[nrow(non_cg_heatmap):1, ]),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Non-CG Contexts Heatmap", axes = FALSE, zlim = c(-1, 1))
  axis(1, at = 1:ncol(non_cg_heatmap), labels = colnames(non_cg_heatmap), las = 2, cex.axis = 0.6)
  axis(2, at = 1:nrow(non_cg_heatmap), labels = rev(rownames(non_cg_heatmap)), las = 2, cex.axis = 0.7)
  
  dev.off()
  
  # 3. Annotation-specific comparison plot
  png("ANNOTATION_SPECIFIC_summary_mto1_correlations.png", 
      width = 1800, height = 1000, res = 120, family = "serif")
  
  layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE))
  
  # Panel 1: Genes across all contexts
  par(mar = c(6, 4, 4, 2))
  genes_data <- all_correlations[all_correlations$annotation == "Genes", ]
  genes_by_context <- split(genes_data$correlation, genes_data$context)
  boxplot(genes_by_context,
          main = "Genes: All Contexts",
          ylab = "Correlation with mto1",
          col = c("lightblue", "lightgreen", "lightcoral"),
          ylim = c(-1, 1))
  abline(h = 0, lty = 2, col = "darkgray")
  
  # Panel 2: Promoters across all contexts
  par(mar = c(6, 4, 4, 2))
  promoters_data <- all_correlations[all_correlations$annotation == "Promoters", ]
  promoters_by_context <- split(promoters_data$correlation, promoters_data$context)
  boxplot(promoters_by_context,
          main = "Promoters: All Contexts", 
          ylab = "Correlation with mto1",
          col = c("lightblue", "lightgreen", "lightcoral"),
          ylim = c(-1, 1))
  abline(h = 0, lty = 2, col = "darkgray")
  
  # Panel 3: TEs across all contexts
  par(mar = c(6, 4, 4, 2))
  tes_data <- all_correlations[all_correlations$annotation == "Transposable_Elements", ]
  tes_by_context <- split(tes_data$correlation, tes_data$context)
  boxplot(tes_by_context,
          main = "Transposable Elements: All Contexts",
          ylab = "Correlation with mto1", 
          col = c("lightblue", "lightgreen", "lightcoral"),
          ylim = c(-1, 1))
  abline(h = 0, lty = 2, col = "darkgray")
  
  # Panel 4: Top genes correlations
  par(mar = c(4, 8, 4, 2))
  top_genes <- genes_data %>%
    group_by(mutant) %>%
    summarise(mean_cor = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(abs(mean_cor))) %>%
    head(10) %>%
    arrange(mean_cor)
  
  barplot(top_genes$mean_cor,
          names.arg = top_genes$mutant,
          main = "Top Gene Correlations",
          xlab = "Mean Correlation",
          col = ifelse(top_genes$mean_cor > 0, "steelblue", "orange"),
          horiz = TRUE, las = 1, xlim = c(-1, 1))
  abline(v = 0, lty = 1, col = "black", lwd = 2)
  
  # Panel 5: Top promoter correlations
  par(mar = c(4, 8, 4, 2))
  top_promoters <- promoters_data %>%
    group_by(mutant) %>%
    summarise(mean_cor = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(abs(mean_cor))) %>%
    head(10) %>%
    arrange(mean_cor)
  
  barplot(top_promoters$mean_cor,
          names.arg = top_promoters$mutant,
          main = "Top Promoter Correlations",
          xlab = "Mean Correlation",
          col = ifelse(top_promoters$mean_cor > 0, "steelblue", "orange"),
          horiz = TRUE, las = 1, xlim = c(-1, 1))
  abline(v = 0, lty = 1, col = "black", lwd = 2)
  
  # Panel 6: Top TE correlations  
  par(mar = c(4, 8, 4, 2))
  top_tes <- tes_data %>%
    group_by(mutant) %>%
    summarise(mean_cor = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(abs(mean_cor))) %>%
    head(10) %>%
    arrange(mean_cor)
  
  barplot(top_tes$mean_cor,
          names.arg = top_tes$mutant,
          main = "Top TE Correlations",
          xlab = "Mean Correlation", 
          col = ifelse(top_tes$mean_cor > 0, "steelblue", "orange"),
          horiz = TRUE, las = 1, xlim = c(-1, 1))
  abline(v = 0, lty = 1, col = "black", lwd = 2)
  
  dev.off()
  
  # Print summary statistics
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("COMPREHENSIVE CORRELATION SUMMARY\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Total combinations analyzed:", nrow(all_correlations), "\n")
  cat("Unique mutants:", length(unique(all_correlations$mutant)), "\n")
  cat("Contexts analyzed:", paste(unique(all_correlations$context), collapse = ", "), "\n")
  cat("Annotations analyzed:", paste(unique(all_correlations$annotation), collapse = ", "), "\n\n")
  
  cat("Summary plots created:\n")
  cat("1. CG_CONTEXT_summary_mto1_correlations.png\n")
  cat("2. NON_CG_CONTEXT_summary_mto1_correlations.png\n") 
  cat("3. ANNOTATION_SPECIFIC_summary_mto1_correlations.png\n")
  
  return(all_correlations)
}

# RUN THE COMPREHENSIVE SUMMARY ANALYSIS
summary_results <- create_comprehensive_summary_plot()
