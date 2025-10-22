DMRs_n_tair_list_plot <- function(tair_ids, DMRs_directory, gene_line_color = NULL, gene_line_color_legend = FALSE, output_prefix = "", output_path = "./", annotations_file = "https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/annotation_files/Methylome.At_annotations.csv.gz", TEs_file = "https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/annotation_files/TAIR10_Transposable_Elements.txt", gene2symbol_file = "https://raw.githubusercontent.com/Yo-yerush/RA_lab_db/main/Arabidopsis/tair2symbol.txt") {
    # -----------------------------------------------
    # run example:
    # tair_ids <- c("AT2G41230", "AT1G36060", "AT3G16770")
    # DMRs_directory <- "PATH/TO/Methylome.At/results/mto1_vs_wt/genome_annotation/"
    # DMRs_n_tair_list_plot(tair_ids, DMRs_directory, gene_line_color = "rainbow", gene_line_color_legend = T)
    # DMRs_Density_legend()
    # -----------------------------------------------

    library(dplyr)
    library(GenomicRanges)
    library(circlize)

    source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/trimm_and_rename_seq.R")
    source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/edit_TE_file.R")

    annotation.gr <- data.table::fread(annotations_file) %>%
        as.data.frame() %>%
        makeGRangesFromDataFrame(keep.extra.columns = T) %>%
        trimm_and_rename()

    TEs.gr <- read.csv(TEs_file, sep = "\t") %>% edit_TE_file()

    chr_amount <- length(seqnames(annotation.gr)@values)

    colors_vec <- if (!is.null(gene_line_color)) {
        gene_line_color
    } else if (length(tair_ids) > 9) {
        rainbow(length(tair_ids))
    } else {
        RColorBrewer::brewer.pal(length(tair_ids), "Set1")
    }

    tair_ids <- toupper(tair_ids)
    tair_list_col <- data.frame(gene_id = tair_ids, col = colors_vec)
    tair_list <- tair_ids %>% paste(collapse = "|")
    tair_list_gr <- annotation.gr[which(annotation.gr$type == "gene")]
    tair_list_gr <- tair_list_gr[grep(tair_list, tair_list_gr$gene_id)]
    start(tair_list_gr) <- (start(tair_list_gr) + end(tair_list_gr)) / 2
    end(tair_list_gr) <- start(tair_list_gr)
    tair_col <- tair_list_col[match(tair_list_gr$gene_id, tair_list_col$gene_id), ]

    cntx_file <- function(context) {
        ############# read DMRs file
        dmrs_file <- read.csv(list.files(paste0(DMRs_directory, "../"), pattern = paste0("DMRs_", context, "_.*\\.csv"), full.names = TRUE))
        dmrs_file <- dmrs_file[, c("seqnames", "start", "end", "log2FC")]
        return(dmrs_file)
    }

    CG_file <- cntx_file("CG")
    CHG_file <- cntx_file("CHG")
    CHH_file <- cntx_file("CHH")

    #####################################
    output_prefix <- ifelse(output_prefix == "", "", paste0(output_prefix, "_"))
    output_prefix <- paste0(output_path, output_prefix)
    ############# the plot #############
    svg(paste0(output_prefix, "tair_list_DMRs_Density.svg"), width = 3.25, height = 3.25, family = "serif")

    circos.par(start.degree = 90)
    circos.genomicInitialize(as.data.frame(annotation.gr)[, 1:3], sector.names = paste0("Chr ", seq(chr_amount)), axis.labels.cex = 0.325, labels.cex = 1.35)

    suppressWarnings({
        circos.genomicDensity(
            list(
                CG_file[CG_file$log2FC > 0, 1:3],
                CG_file[CG_file$log2FC < 0, 1:3]
            ),
            bg.col = "#fafcff", bg.border = NA, count_by = "number",
            col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
        )

        circos.genomicDensity(
            list(
                CHG_file[CHG_file$log2FC > 0, 1:3],
                CHG_file[CHG_file$log2FC < 0, 1:3]
            ),
            bg.col = "#fafcff", bg.border = NA, count_by = "number",
            col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
        )

        circos.genomicDensity(
            list(
                CHH_file[CHH_file$log2FC > 0, 1:3],
                CHH_file[CHH_file$log2FC < 0, 1:3]
            ),
            bg.col = "#fafcff", bg.border = NA, count_by = "number",
            col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
        )

        circos.genomicTrackPlotRegion(as.data.frame(tair_list_gr)[, 1:3], ylim = c(0, 1), bg.col = "#fafcfc", bg.border = NA, panel.fun = function(region, value, ...) {
            region_colors <- tair_col$col[which(as.data.frame(tair_list_gr)[, 1] == get.current.chromosome())]
            circos.genomicLines(region, 1, type = "h", col = region_colors, lwd = 0.8)
        }, track.height = 0.06, track.margin = c(0, 0))

        circos.genomicDensity(as.data.frame(TEs.gr)[, 1:3],
            bg.col = "#fafcff", bg.border = NA, count_by = "number",
            col = "#fcba0320", border = T, track.height = 0.165, track.margin = c(0, 0)
        )

        circos.clear()
        dev.off()
    })

    if (gene_line_color_legend == T) {
        # tair_list_col_ordered <- tair_list_col[match(tair_list_gr$gene_id, tair_list_col$gene_id), ]
        tair_list_col_ordered <- merge(data.frame(gene_id=tair_list_gr$gene_id), tair_list_col, by = "gene_id", all.x = T)

        gene2symbol <- read.csv(gene2symbol_file, sep = "\t", header = F)
        names(gene2symbol) <- c("gene_id", "Symbol")
        symbol_match <- match(tair_list_col_ordered$gene_id, gene2symbol[, 1])
        tair_list_col_ordered$symbol <- gene2symbol[symbol_match, 2]

        plot_height <- max(2.5, ceiling(nrow(tair_list_col_ordered) / 4) * 1.25)
        top_margin <- 2.5

        pdf(paste0(output_prefix, "tair_color_legend.pdf"), width = 3, height = plot_height, family = "serif")
        par(mar = c(0.1, 0.1, top_margin, 0.1))
        plot(1,
            type = "n", xlim = c(0, 10), ylim = c(0, nrow(tair_list_col_ordered)),
            xlab = "", ylab = "", axes = FALSE
        )
        for (i in 1:nrow(tair_list_col_ordered)) {
            y_pos <- nrow(tair_list_col_ordered) - i + 1

            text(0.25, y_pos, tair_list_col_ordered$gene_id[i],
                adj = c(0, 0.5), cex = 0.9
            )

            if ("symbol" %in% colnames(tair_list_col_ordered) && tair_list_col_ordered$symbol[i] != "") {
                text(7.15, y_pos, tair_list_col_ordered$symbol[i],
                    adj = c(0, 0.5), cex = 0.9, col = "black"
                )
            }

            rect(3.5, y_pos - 0.2, 6.75, y_pos + 0.2,
                col = tair_list_col_ordered$col[i], border = "black"
            )
        }
        title("TAIR IDs and Colors\n*ordered by genome position", cex.main = 1.1)
        dev.off()
    }

    # output the gene list data frame as vector
    as.data.frame(tair_list_gr) %>%
    merge(., gene2symbol, by = "gene_id", all.x = T) %>%
    merge(., tair_list_col, by = "gene_id", all.x = T) %>%
    dplyr::relocate(gene_id, Symbol, col, .after = strand) %>%
    select(-end, -width, -gene_model_type)
}
############################################################

# legends
DMRs_Density_legend <- function() {
    rndm_dis <- function() {
        counts_norm_hyper <- hist(rnorm(3000, mean = 0, sd = 1), plot = F, breaks = 30)$counts
        counts_norm_hypo <- hist(rnorm(3000, mean = 0, sd = 1), plot = F, breaks = 30)$counts
        rndm_norm_hyper <- c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0, 10))
        rndm_norm_hypo <- c(rep(0, 10), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
        # rndm_norm_TE = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
        # rndm_norm_TE[rndm_norm_TE < 0.05] = rnorm(length(rndm_norm_TE[rndm_norm_TE < 0.05]), mean = 0.05, sd = 0.02)
        # rndm_norm_TE = rndm_norm_TE[!rndm_norm_TE == 0]
        return(list(hyper = rndm_norm_hyper, hypo = rndm_norm_hypo))
    }
    hyper.1 <- rndm_dis()$hyper
    hyper.2 <- rndm_dis()$hyper
    hyper.3 <- rndm_dis()$hyper
    hypo.1 <- rndm_dis()$hypo
    hypo.2 <- rndm_dis()$hypo
    hypo.3 <- rndm_dis()$hypo

    rndm_norm_TE <- c(0, rndm_dis()$hypo[-c(1:10)], 0)

    rndm <- runif(12, min = 0, max = 1)
    legend_titles_circular <- c("DMRs", "  Group sig. genes", "Methylation values", "Expression values", "  Overlapping TEs", "Genome annotation")
    tracks_total_size <- 16
    legend_tracks_pos <- as.character((tracks_total_size / 4) + 1)
    text_tracks_pos <- as.character((tracks_total_size / 4))

    colfunc <- colorRampPalette(c("red", "white", "blue"))
    colfunc_vec_1 <- colfunc(16)[-c((round(16 / 2) - 1):(16 / 2), ((16 / 2) + 1):((16 / 2) + 1 + 1))]
    colfunc_vec_2 <- colfunc(100)[-c((round(100 / 2) - 5):(100 / 2) - 1, ((100 / 2) + 1):((100 / 2) + 1 + 5))]


    svg("legends.svg", width = 2.5, height = 3.35, family = "serif")
    par(mar = c(0, 0, 0, 0))
    circos.par("track.height" = 0.6, "canvas.xlim" = c(-0.2, 0.3), "canvas.ylim" = c(-0.25, 1), "gap.degree" = 0, "clock.wise" = FALSE)
    circos.initialize(factors = as.character(1:tracks_total_size), xlim = c(0, 1))

    circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175)
    circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA, bg.col = "#fafcfc", bg.lwd = 0.2)
    circos.lines(seq(0, 1, 1 / (length(hyper.1) - 1)), hyper.1, area = T, border = T, col = "#FF000080")
    circos.lines(seq(0, 1, 1 / (length(hypo.1) - 1)), hypo.1, area = T, border = T, col = "#304ed180")
    circos.lines(seq(0, 1, 1 / (length(hyper.1) - 1)), hyper.1, border = T)
    circos.text(text_tracks_pos, x = 0.1, y = 0.5, labels = "CG context", facing = "downward", ad = c(0, 0.5))

    circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175)
    circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA, bg.col = "#fafcfc", bg.lwd = 0.2)
    circos.lines(seq(0, 1, 1 / (length(hyper.2) - 1)), hyper.2, area = T, border = T, col = "#FF000080")
    circos.lines(seq(0, 1, 1 / (length(hypo.2) - 1)), hypo.2, area = T, border = T, col = "#304ed180")
    circos.lines(seq(0, 1, 1 / (length(hyper.2) - 1)), hyper.2, border = T)
    circos.text(text_tracks_pos, x = 0.14, y = 0.5, labels = "CHG context", facing = "downward", adj = c(0, 0.5))

    circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175)
    circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA, bg.col = "#fafcfc", bg.lwd = 0.2)
    circos.lines(seq(0, 1, 1 / (length(hyper.3) - 1)), hyper.3, area = T, border = T, col = "#FF000080")
    circos.lines(seq(0, 1, 1 / (length(hypo.3) - 1)), hypo.3, area = T, border = T, col = "#304ed180")
    circos.lines(seq(0, 1, 1 / (length(hyper.3) - 1)), hyper.3, border = T)
    circos.text(text_tracks_pos, x = 0.22, y = 0.5, labels = "CHH context", facing = "downward", adj = c(0, 0.5))

    circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.15)
    circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA, bg.col = "#fafcfc", bg.lwd = 0.2)
    circos.lines(seq(0, 1, 1 / (length(rndm_norm_TE) - 1)), rndm_norm_TE, area = T, border = T, col = "#fcba0330")
    circos.text(text_tracks_pos, x = 0.37, y = 0.5, labels = "TEs", facing = "downward", adj = c(0, 0.5))

    circos.clear()

    # DMRs
    legend(-0.12, 0.125,
        legend = c(
            substitute(paste(bold("Hyper"), "-DMRs")),
            substitute(paste(bold("Hypo"), "-DMRs")),
            "Shared"
        ),
        fill = c("#FF000095", "#304ed195", "#5e0d3d99"),
        bty = "n"
    )

    dev.off()
}
