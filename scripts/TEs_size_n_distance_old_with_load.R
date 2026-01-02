suppressMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(DMRcaller)
    library(org.At.tair.db)
    library(GenomicFeatures)
    library(plyranges)
    library(parallel)
    library(data.table)
})

##################################################################################

source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/trimm_and_rename_seq.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/load_replicates.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/edit_TE_file.R")

##################################################################################

message("create TE objects into env")
TE_file_path <- "https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/annotation_files/TAIR10_Transposable_Elements.txt"
TE_df <- read.csv(TE_file_path, sep = "\t")
te_width <- data.frame(
    te_id = TE_df$Transposon_Name,
    width = (TE_df$Transposon_max_End - TE_df$Transposon_min_Start),
    superfamily = (TE_df$Transposon_Super_Family)
)
long_tes <- te_width[te_width$width >= 4000, 1]
short_tes <- te_width[te_width$width <= 500, 1]

TE_gr <- edit_TE_file(TE_df)

rm(list = c("TE_df", "TE_file_path"))

##################################################################################

TE_delta_meth <- function(var_table = NULL, pooled_list = NULL, context = NULL, n.cores = 6) {

    if (!is.null(pooled_list)) {
        # use pooled and trimmed object
        meth_ctrl <- pooled_list[[1]]
        meth_trnt <- pooled_list[[2]]

    } else if (!is.null(var_table)) {
        # load replicates raw 'CX' data
        vars_vector <- unique(var_table[, 1])
        var1_name <- vars_vector[1]
        var2_name <- vars_vector[2]
        var1_path <- var_table[grep(var1_name, var_table[, 1]), 2]
        var2_path <- var_table[grep(var2_name, var_table[, 1]), 2]
        var_args <- list(
            list(path = var1_path, name = var1_name),
            list(path = var2_path, name = var2_name)
        )
        load_vars <- mclapply(var_args, function(x) load_replicates(x$path, ifelse(n.cores > 1, n.cores / 2, 1), x$name, T, "CX_report"), mc.cores = ifelse(n.cores > 1, 2, 1))
        meth_ctrl <- trimm_and_rename(load_vars[[1]])
        meth_trnt <- trimm_and_rename(load_vars[[2]])
    }

    # Calculate average methylation for both samples
    meth_delta <- meth_ctrl[, 1]
    meth_delta$Proportion <- (meth_trnt$readsM / meth_trnt$readsN) - (meth_ctrl$readsM / meth_ctrl$readsN)
    meth_delta$Proportion[is.nan(meth_delta$Proportion)] <- 0

    # plots
    te_meth_delta <- list()
    for (cntx in c("CG", "CHG", "CHH")) {
        # delta
        te_meth_delta[[cntx]] <- calculate_te_methylation(meth_delta, TE_gr, cntx, T)
        te_meth_delta[[cntx]]$sample <- "delta"
    }

    return(te_meth_delta)
}

##################################################################################

calculate_te_methylation <- function(meth_data, TE_gr, context, is.delta = F) {
    meth_data <- meth_data[meth_data$context == context]
    if (!is.delta) {
        meth_data$Proportion <- meth_data$readsM / meth_data$readsN
    }
    # Find overlaps between methylation data and TEs
    overlaps <- findOverlaps(TE_gr, meth_data)
    overlaps_df <- data.frame(
        te_id = queryHits(overlaps),
        meth_idx = subjectHits(overlaps)
    )

    # Get TE IDs and methylation levels for overlapping regions
    overlaps_df$te_id <- TE_gr$gene_id[overlaps_df$te_id]
    overlaps_df$meth_level <- meth_data$Proportion[overlaps_df$meth_idx]

    # Calculate average methylation level for each TE
    te_avg_meth <- overlaps_df %>%
        group_by(te_id) %>%
        summarise(avg_meth = mean(meth_level, na.rm = TRUE), .groups = "drop")
    te_avg_meth$context <- context

    filtered_te_width <- te_width %>% filter(width < 20000)

    as.data.frame(te_avg_meth) %>% merge(., filtered_te_width, by = "te_id")
}

##################################################################################

te_size_plot <- function(x, cntx, line_col = "red4", point_col = "gray20") {
    ggplot(x[[cntx]], aes(x = width, y = avg_meth, color = sample)) +
        geom_vline(xintercept = 4000, linetype = "dashed", color = "gray60") +
        geom_point(alpha = 0.6, size = 0.3, shape = 20) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        geom_smooth(
            method = lm, formula = y ~ splines::bs(x, 3),
            se = TRUE,
            linewidth = 0.5,
            color = line_col
        ) +
        scale_color_manual(values = c("delta" = point_col)) +
        labs(
            title = paste(cntx, "context"),
            x = "TE size (Kbp)",
            y = paste0("Δ ", "Methylation")
        ) +
        theme_classic() +
        theme(
            text = element_text(family = "serif"),
            legend.position = "none",
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            plot.title = element_text(hjust = 0.5, size = 10),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 9)
        ) +
        scale_x_continuous(
            breaks = c(65, 5000, 10000, 15000, 19935),
            labels = seq(0, 20, by = 5),
            limits = c(0, 20000), expand = c(0, 0)
        ) +
        scale_y_continuous(limits = c(-0.3, 0.3), expand = c(0, 0), breaks = c(-0.296, 0, 0.296), labels = c(-0.3, 0, 0.3))
}

##################################################################################

distance_from_centromer <- function(TE_meth_delta_list, TE_gr = TE_gr, window_size = 1e6, lines_col = c("#3d53b4", "#3b8f3e", "#bb4949"), y_max = NULL, y_min = NULL, y_breaks = NULL) {
    # centromers positions
    cen_pos <- c(14.845, 3.44, 13.855, 3.13, 11.795) * 1e6
    te_distance <- data.frame(
        chr = as.character(seqnames(TE_gr)),
        pos = (as.numeric(start(TE_gr)) + as.numeric(end(TE_gr))) / 2,
        te_id = TE_gr$gene_id,
        centromere = NA
    )
    for (cen_i in 1:5) {
        te_distance$centromere[te_distance$chr == paste0("Chr", cen_i)] <- cen_pos[cen_i]
    }
    te_distance$distance <- abs(te_distance$centromere - te_distance$pos)
    all_cx_dis <- rbind(
        TE_meth_delta_list[["CG"]],
        TE_meth_delta_list[["CHG"]],
        TE_meth_delta_list[["CHH"]]
    )
    keep_col <- grep("te_id|distance", names(te_distance))
    te_distance_merged_out <- merge(te_distance[, keep_col], all_cx_dis, by = "te_id")
    te_distance_merged <- te_distance_merged_out
    te_distance_merged$distance <- te_distance_merged$distance / window_size

    ########

    te_distance_cntx <- rbind(
        te_distance_merged %>%
            filter(context == "CG") %>%
            mutate(window = floor(distance)) %>% # floor(distance / 1e6) * 1e6) %>%
            group_by(window, context) %>%
            summarise(
                avg_meth = mean(avg_meth, na.rm = TRUE),
                distance = mean(distance, na.rm = TRUE),
                .groups = "drop"
            ),
        te_distance_merged %>%
            filter(context == "CHG") %>%
            mutate(window = floor(distance)) %>% # floor(distance / 1e6) * 1e6) %>%
            group_by(window, context) %>%
            summarise(
                avg_meth = mean(avg_meth, na.rm = TRUE),
                distance = mean(distance, na.rm = TRUE),
                .groups = "drop"
            ),
        te_distance_merged %>%
            filter(context == "CHH") %>%
            mutate(window = floor(distance)) %>% # floor(distance / 1e6) * 1e6) %>%
            group_by(window, context) %>%
            summarise(
                avg_meth = mean(avg_meth, na.rm = TRUE),
                distance = mean(distance, na.rm = TRUE),
                .groups = "drop"
            )
    ) %>%
        as.data.frame() %>%
        filter(avg_meth != max(avg_meth, na.rm = TRUE))

    ########

    y_axis_limits <- if (!is.null(y_max) | !is.null(y_min)) {
        if (!is.null(y_breaks)) {
            y_breaks <- c(y_min, 0, y_max)
        }
        scale_y_continuous(
            limits = c(y_min, y_max),
            breaks = y_breaks
        )
    } else {
        geom_blank()
    }

    ########

    te_distance_plot <- ggplot(data = te_distance_cntx, aes(x = distance, y = avg_meth, color = context, group = context)) +
        geom_line(linewidth = 0.85) +
        theme_bw() + # theme_classic() +
        labs(
            title = " ",
            x = "Distance from centromer (Mbp)",
            y = "Δ Methylation"
        ) +
        theme(
            text = element_text(family = "serif"),
            legend.position = "none",
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            plot.title = element_text(hjust = 0.5, size = 10),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 9)
        ) +
        scale_x_continuous(
            limits = c(0, 15), # max(-te_distance_cntx$distance)),
            breaks = c(0.06, 5, 10, 14.94),
            labels = seq(0, 15, by = 5),
            expand = c(0, 0)
        ) +
        y_axis_limits +
        annotate("text",
            x = 12.75, # 12.25,
            y = max(te_distance_cntx$avg_meth) * 0.98,
            label = c("CG", "\nCHG", "\n\nCHH"),
            hjust = 0, vjust = 0.75, size = 3.25,
            color = lines_col, fontface = "bold", family = "serif"
        )

    return(list(
        df = select(te_distance_merged_out, -sample),
        plot = te_distance_plot
    ))
}
