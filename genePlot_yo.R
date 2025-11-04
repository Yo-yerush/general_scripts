# Yo - 080725
#
# -------------------------------------------------------------

# DMRcaller packge based functions
#
# Catoni, Marco, Tsang, MF J, Greco, P A, Zabet, Radu N (2018). “DMRcaller: a versatile R/Bioconductor package for detection and visualization of differentially methylated regions in CpG and non-CpG contexts.” Nucleic Acids Research. doi:10.1093/nar/gky602.
#
# https://www.bioconductor.org/packages/release/bioc/html/DMRcaller.html

# -------------------------------------------------------------
#
# script to plot methylation level over gene body and promoter
# insters TAIR ID(s), pooled CX_reports and variable names
#
# 'methylome_at_annotations' in default use TAIR10 annotations (both genes and TEs)

# -------------------------------------------------------------

run_DMRs_genePlots <- function(tair_id, ctrl_gr, mut_gr, ctrl_name, mut_name, output_path = ".", DMP_calling = TRUE, DMP_min_diff = c(40, 20, 10), DMP_qValue = 0.05, DMP_min_cov = 4, DMPs_table = FALSE, DMPs_pie = FALSE, nbp_upstream = 2000, nbp_downstream = 200, genePlot_legend = TRUE) {
    suppressMessages({
        library(DMRcaller)
        library(GenomicFeatures)
        library(dplyr)
        library(methylKit)
    })
    exp_condition <- paste0(mut_name, "_vs_", ctrl_name)

    ###################################

    cat("\rupload annotation files...     ")
    annotations_url <- "https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/annotation_files/Methylome.At_annotations.csv.gz"
    description_url <- "https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/annotation_files/Methylome.At_description_file.csv.gz"
    transposons_url <- "https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/annotation_files/TAIR10_Transposable_Elements.txt"

    ann_file <- read.csv(
        gzcon(url(annotations_url, open = "rb"), text = TRUE),
        encoding = "UTF-8"
    ) %>%
        dplyr::select(-width)

    desc_file <- read.csv(
        gzcon(url(description_url, open = "rb"), text = TRUE),
        encoding = "UTF-8"
    ) %>%
        dplyr::select(gene_id, Symbol)

    TE_file <- read.csv(
        gzcon(url(transposons_url, open = "rb"), text = TRUE),
        encoding = "UTF-8",
        sep = "\t"
    ) %>%
        mutate(seqnames = NA) %>%
        mutate(type = "transposable_element", gene_model_type = "transposable_element") %>%
        dplyr::select(seqnames, Transposon_min_Start, Transposon_max_End, orientation_is_5prime, type, Transposon_Name, gene_model_type) %>%
        dplyr::rename(gene_id = Transposon_Name)
    for (i in 1:5) {
        TE_file$seqnames[grep(paste0("AT", i, "TE"), TE_file$gene_id)] <- paste0("Chr", i)
    }
    TE_file$orientation_is_5prime <- gsub("true", "+", TE_file$orientation_is_5prime)
    TE_file$orientation_is_5prime <- gsub("false", "-", TE_file$orientation_is_5prime)
    names(TE_file)[1:4] <- c("seqnames", "start", "end", "strand")

    gff3_file <- rbind(ann_file, TE_file) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    cat("\rupload annotation files... done")
    cat("\n\n")

    ###################################

    # get gene position from the gff3 file
    gene_gr <- as.data.frame(gff3_file) %>%
        filter(type == "gene") %>%
        filter(gene_id == tair_id) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    if (as.character(strand(gene_gr)) == "+") {
        start(gene_gr) <- start(gene_gr) - nbp_upstream
        end(gene_gr) <- end(gene_gr) + nbp_downstream
    } else {
        start(gene_gr) <- start(gene_gr) - nbp_downstream
        end(gene_gr) <- end(gene_gr) + nbp_upstream
    }

    chr_name <- as.character(gene_gr@seqnames@values[1])
    start_pos <- start(gene_gr)
    end_pos <- end(gene_gr)

    # Filter var_pool by the positions of chr_name, start_pos, and end_pos
    filtered_ctrl <- ctrl_gr[seqnames(ctrl_gr) == chr_name & start(ctrl_gr) >= start_pos & end(ctrl_gr) <= end_pos]
    filtered_mut <- mut_gr[seqnames(mut_gr) == chr_name & start(mut_gr) >= start_pos & end(mut_gr) <= end_pos]

    joint_gr <- .joinMethylationData(ctrl_gr, mut_gr)
    condition <- c(ctrl_name, mut_name) ### fix

    ###################################

    if (DMP_calling) {
        ### DMPs - using methylkit
        dmps_CG <- methylkit_dmps_joint_as_gr(joint_gr, condition, cntx = "CG", diff_main = DMP_min_diff[1], q_main = DMP_qValue, min_cov = DMP_min_cov)
        dmps_CHG <- methylkit_dmps_joint_as_gr(joint_gr, condition, cntx = "CHG", diff_main = DMP_min_diff[2], q_main = DMP_qValue, min_cov = DMP_min_cov)
        dmps_CHH <- methylkit_dmps_joint_as_gr(joint_gr, condition, cntx = "CHH", diff_main = DMP_min_diff[3], q_main = DMP_qValue, min_cov = DMP_min_cov)
        DMPs_calc_mk <- c(dmps_CG, dmps_CHG, dmps_CHH)
        DMPs_calc_mk <- DMPs_calc_mk[!((DMPs_calc_mk$sumReadsM1 + DMPs_calc_mk$sumReadsM2) == 0)]
        DMPs_calc_mk <- DMPs_calc_mk[!((DMPs_calc_mk$sumReadsM1 + DMPs_calc_mk$sumReadsM2) == 2)]
        # DMPs_calc_mk <- DMPs_calc_mk[DMPs_calc_mk$pValue <= 0.05]
        DMPs_calc_mk <- DMPs_calc_mk[order(start(DMPs_calc_mk))]

        cg <- DMPs_calc_mk[DMPs_calc_mk$context == "CG"]
        chg <- DMPs_calc_mk[DMPs_calc_mk$context == "CHG"]
        chh <- DMPs_calc_mk[DMPs_calc_mk$context == "CHH"]

        DMRsList <- list("CG" = cg, "CHG" = chg, "CHH" = chh)

        if (create_legend) {
            legend_genePlot_DMPs(output_path)
        }
    } else if (DMRs_list) {

    } else {
        DMPs_calc_mk <- NULL
    }
    ###################################

    # DMRs/DMPs
    if (!is.null(MPs_calc_mk)) {

    } else {
        DMRsList <- NULL
    }

    if (DMPs_table) {
        write.csv(DMPs_calc_mk, paste0(output_path, "/", tair_id, "_DMPs_", exp_condition, ".csv"), row.names = F)
    }

    if (DMPs_pie) {
        gainORloss(DMPs_calc_mk, TRUE, paste0(exp_condition, "/using_methylKit"))
    }

    ###################################


    genePlot_fun(tair_id, ctrl_gr, mut_gr, ctrl_name, mut_name, DMRsList = DMRsList, output_path = output_path, create_legend = genePlot_legend)

    #############
}

################################################################################

genePlot_fun <- function(tair_id, var1_pool, var2_pool, var1_name, var2_name, DMRsList = NULL, output_path = ".", create_legend = FALSE, n_bp_upstream = 2000, n_bp_downstream = 1000) {
    ###################################

    id2symbol <- filter(desc_file, gene_id == tair_id)

    main_title <- ifelse(
        is.na(id2symbol$Symbol),
        id2symbol$gene_id,
        paste0(id2symbol$gene_id, " (", id2symbol$Symbol, ")")
    )
    file_suffix <- gsub("\\(|\\)", "", main_title)
    file_suffix <- gsub(" ", "_", file_suffix)
    file_suffix <- gsub("::", "_", file_suffix)

    ###################################

    svg(paste0(output_path, "/genePlot_", file_suffix, ".svg"), width = 4.75, height = 2.75, family = "serif")
    plotLocalMethylationProfile_yo(
        filtered_var1_pool,
        filtered_var2_pool,
        gene_gr,
        DMRsList,
        conditionsNames = c(var1_name, var2_name),
        gff3_file,
        context = c("CG", "CHG", "CHH"),
        main = main_title,
        col = c(
            "#404040", "#bf6828", "#000000", "#afad43",
            "#009E73", "#0072B2", "#CC79A7", "#D55E00", "#999999"
        )
    )
    dev.off()
}

##########################################
### functions
# rename_seq <- function(gr_obj) {
#     gr_obj <- renameSeqlevels(gr_obj, gsub("NC_003070.9", "Chr1", seqlevels(gr_obj)))
#     gr_obj <- renameSeqlevels(gr_obj, gsub("NC_003071.7", "Chr2", seqlevels(gr_obj)))
#     gr_obj <- renameSeqlevels(gr_obj, gsub("NC_003074.8", "Chr3", seqlevels(gr_obj)))
#     gr_obj <- renameSeqlevels(gr_obj, gsub("NC_003075.7", "Chr4", seqlevels(gr_obj)))
#     gr_obj <- renameSeqlevels(gr_obj, gsub("NC_003076.8", "Chr5", seqlevels(gr_obj)))
#     return(gr_obj)
# }

.joinMethylationData <- function(cx1, cx2) {
    overlaps <- findOverlaps(cx1, cx2)
    indexes <- which(!duplicated(queryHits(overlaps)))
    methylData <- GRanges(
        seqnames = seqnames(cx1[queryHits(overlaps)[indexes]]),
        ranges = ranges(cx1[queryHits(overlaps)[indexes]]),
        strand = strand(cx1[queryHits(overlaps)[indexes]]),
        context = cx1$context[queryHits(overlaps)[indexes]],
        trinucleotide_context = cx1$trinucleotide_context[queryHits(overlaps)[indexes]],
        readsM1 = cx1$readsM[queryHits(overlaps)[indexes]],
        readsN1 = cx1$readsN[queryHits(overlaps)[indexes]],
        readsM2 = cx2$readsM[subjectHits(overlaps)[indexes]],
        readsN2 = cx2$readsN[subjectHits(overlaps)[indexes]]
    )
    return(methylData)
}

.plotGeneticElements <- function(gff, region, col) {
    seqname <- seqnames(region)
    minPos <- start(region)
    maxPos <- end(region)


    # Select the genes that lie in the region of interest.
    gff <- gff[queryHits(findOverlaps(gff, region))]
    # Chop off the ends of anything sticking out...
    start(gff) <- pmax(start(gff), minPos)
    end(gff) <- pmin(end(gff), maxPos)

    genes <- gff[gff$type == "gene"]
    genesPos <- genes[strand(genes) == "+" | strand(genes) == "*"]
    genesNeg <- genes[strand(genes) == "-" | strand(genes) == "*"]
    exons <- gff[gff$type == "exon"]
    exons <- exons[overlapsAny(exons, genes)]
    exonsPos <- exons[strand(exons) == "+" | strand(exons) == "*"]
    exonsNeg <- exons[strand(exons) == "-" | strand(exons) == "*"]

    transposons <- gff[gff$type == "transposable_element"]
    transposonsPos <- transposons[strand(transposons) == "+" | strand(transposons) == "*"]
    transposonsNeg <- transposons[strand(transposons) == "-" | strand(transposons) == "*"]

    negativeStrandPosition <- -0.175
    positiveStrandPosition <- -0.075


    # text(maxPos + (maxPos-minPos)/100, positiveStrandPosition, 'Sense', cex = 0.5);
    # text(maxPos + (maxPos-minPos)/100, negativeStrandPosition, 'Anti', cex = 0.5);
    text(maxPos + (maxPos - minPos) / 100, positiveStrandPosition, "+", font = 2)
    text(maxPos + (maxPos - minPos) / 100, negativeStrandPosition, "-", font = 2)
    lines(c(minPos, maxPos), c(-0.14, -0.14), lty = 1, lwd = 0.75, col = "black")

    if (length(genesPos) > 0) {
        segments(start(genesPos), positiveStrandPosition, end(genesPos), positiveStrandPosition)
        text(start(genesPos), -0.115, genesPos$ID, pos = 4, cex = 0.5)
    }
    if (length(genesNeg) > 0) {
        segments(start(genesNeg), -0.175, end(genesNeg), negativeStrandPosition)
        text(start(genesNeg), -0.23, genesNeg$ID, pos = 4, cex = 0.5)
    }
    if (length(exonsPos) > 0) {
        rect(start(exonsPos), -0.05, end(exonsPos), -0.09, col = col[1], border = NA)
    }
    if (length(exonsNeg) > 0) {
        rect(start(exonsNeg), -0.16, end(exonsNeg), -0.2, col = col[1], border = NA)
    }
    if (length(transposonsPos) > 0) {
        rect(start(transposonsPos), -0.05, end(transposonsPos), -0.09, col = col[2], border = col[2], density = 30, angle = 30)
        text(start(transposonsPos), -0.115, transposonsPos$ID, pos = 4, cex = 0.5)
    }
    if (length(transposonsNeg) > 0) {
        rect(start(transposonsNeg), -0.16, end(transposonsNeg), -0.2, col = col[2], border = col[2], density = 30, angle = 30)
        text(start(transposonsNeg), -0.23, transposonsNeg$ID, pos = 4, cex = 0.5)
    }
}

plotLocalMethylationProfile_yo <- function(
    methylationData, region, DMRs = NULL,
    conditionsNames = NULL, gff = NULL, context = "CG",
    labels = NULL, col = NULL, main = "", plotMeanLines = TRUE,
    plotPoints = TRUE) {
    numberOfConditions <- 2
    if (is.null(conditionsNames) | length(conditionsNames) <
        numberOfConditions | !all(is.character(conditionsNames))) {
        conditionsNames <- paste("condition ", (1:numberOfConditions),
            sep = ""
        )
    }

    numberOfDMRs <- 0
    if (!is.null(DMRs)) {
        numberOfDMRs <- length(DMRs)
    }
    if (!is.null(labels) & (length(labels) < 1 | !is.character(labels))) {
        labels <- LETTERS[1:length(labels)]
    }

    cond1Color <- col[1]
    cond2Color <- col[2]
    geneColor <- col[3]
    TEColor <- col[4]
    DMRsColor <- col[5:length(col)]

    chr_name <- as.character(region@seqnames@values[1])
    seqname <- seqnames(region)
    minPos <- start(region)
    maxPos <- end(region)

    # hits <- findOverlaps(methylationData, region)
    # localMethylationData <- methylationData[queryHits(hits)]
    # contextMethylationData <- localMethylationData[localMethylationData$context %in% context]
    contextMethylationData <- methylationData[methylationData$context %in% context]
    ramp1 <- colorRampPalette(c("white", cond1Color))
    colramp1 <- ramp1(100)
    ramp2 <- colorRampPalette(c("white", cond2Color))
    colramp2 <- ramp2(100)
    proportion1 <- rep(0, length(contextMethylationData))
    index <- which(contextMethylationData$readsM1 >= 0 & contextMethylationData$readsN1 >
        0)
    proportion1[index] <- contextMethylationData$readsM1[index] / contextMethylationData$readsN1[index]
    proportion2 <- rep(0, length(contextMethylationData))
    index <- which(contextMethylationData$readsM2 >= 0 & contextMethylationData$readsN2 >
        0)
    proportion2[index] <- contextMethylationData$readsM2[index] / contextMethylationData$readsN2[index]
    maxColor <- max(c(
        contextMethylationData$readsN1[!is.na(contextMethylationData$readsN1)],
        contextMethylationData$readsN2[!is.na(contextMethylationData$readsN2)]
    ))
    par(mar = c(4, 4, 0, 1) + 0.1)
    plot(start(contextMethylationData), proportion1,
        col = colramp1[round(99 *
            log(contextMethylationData$readsN1) / log(maxColor)) +
            1], pch = 16, cex = 0.6, xlim = c(minPos, maxPos), ylim = c(
            -0.2,
            1.2 + numberOfDMRs * 0.1
        ),
        xlab = paste0("Genomic coordinate (", gsub("Chr", "chromosome ", chr_name), ")"),
        ylab = "Methylation level",
        yaxt = "n", main = NULL, type = "n"
    )
    axis(2, c(0, 0.5, 1))
    if (plotPoints) {
        points(start(contextMethylationData), proportion1, col = colramp1[round(99 *
            log(contextMethylationData$readsN1) / log(maxColor)) +
            1], pch = 16, cex = 0.6)
        points(start(contextMethylationData), proportion2, col = colramp2[round(99 *
            log(contextMethylationData$readsN2) / log(maxColor)) +
            1], pch = 16, cex = 0.6)
        if (!plotMeanLines) {
            # legend("topright", bty = "n", col = c(cond1Color,
            #    cond2Color), legend = conditionsNames, pch = c(16,
            #    16), horiz = TRUE)
        }
    }

    if (numberOfDMRs > 0) {
        for (i in 1:numberOfDMRs) {
            range <- DMRs[[i]][queryHits(findOverlaps(
                DMRs[[i]],
                region
            ))]
            if (length(range) > 0) {
                start(range) <- pmax(start(range), start(region))
                end(range) <- pmin(end(range), end(region))
                if (length(which(range$regionType == "gain")) >
                    0) {
                    rect(start(range)[range$regionType == "gain"],
                        1.075 + (length(DMRs) - i) * 0.1, end(range)[range$regionType ==
                            "gain"], 1.125 + (length(DMRs) - i) * 0.1,
                        col = DMRsColor[i], border = DMRsColor[i]
                    )
                }
                if (length(which(range$regionType == "loss")) >
                    0) {
                    rect(start(range)[range$regionType == "loss"],
                        1.075 + (length(DMRs) - i) * 0.1, end(range)[range$regionType ==
                            "loss"], 1.125 + (length(DMRs) - i) * 0.1,
                        col = DMRsColor[i], border = DMRsColor[i],
                        density = 30
                    )
                }
                if (length(range) > 0 & (is.null(range$regionType))) {
                    rect(start(range), (1.075 + (length(DMRs) -
                        i) * 0.1), end(range), (1.125 + (length(DMRs) -
                        i) * 0.1), col = DMRsColor[i], border = DMRsColor[i])
                }
            }
            ### DMRs in-gragh context labels
            # text(par("usr")[1], 1.1 + (length(DMRs) - i) * 0.1,
            #  pos = 4, names(DMRs)[i], cex = 0.7
            # )
        }
    }

    if (!is.null(gff)) {
        if (is(gff, "GRanges")) {
            .plotGeneticElements(gff, region, c(geneColor, TEColor))
        }
    }
    if (!is.null(labels)) {
        if (length(labels) >= 1) {
            mtext(labels[1], line = 0.7, adj = 0, cex = 1)
        }
    }

    ########################################
    #### add main title
    # text((par("usr")[1] + par("usr")[2]) / 2, 1.45, pos = NULL, main, cex = 0.75, font = 2)
    text(par("usr")[1], 1.45, pos = 4, main, cex = 0.75, font = 2)

    #### DMRs legend
    # text(par("usr")[1], 1.45, pos = 4, "DMRs", cex = 0.7, font = 2)

    #### add legend
    leg_pos1 <- end(region) - (end(region) - start(region)) / 10
    leg_pos2 <- end(region) + (end(region) - start(region)) / 75
    text(leg_pos1, 1.45, pos = 2, conditionsNames[1], cex = 0.75, font = 2, col = cond1Color)
    text(leg_pos2, 1.45, pos = 2, conditionsNames[2], cex = 0.75, font = 2, col = cond2Color)
    # text(par("usr")[2], 1.45, pos = 2, conditionsNames[1], cex = 0.75, font = 2, col = cond1Color)
    # text(par("usr")[2], 1.325, pos = 2, conditionsNames[2], cex = 0.75, font = 2, col = cond2Color)

    #### add line
    # lines(c(
    #  par("usr")[1], par("usr")[2]
    #  #par("usr")[1] - (par("usr")[2] - par("usr")[1]) / 100,
    #  #par("usr")[2] + (par("usr")[2] - par("usr")[1]) / 100
    # ),
    # y = c(1,1),
    #  lty = 1,
    #  lwd = 0.75,
    #  col = "gray75"
    # )

    ########################################

    invisible(NULL)
}

###

#### legend
legend_genePlot <- function(x_path) {
    x.cord <- 0.25
    y.cord <- 0.9

    svg(paste0(x_path, "/genePlot_legend.svg"), width = 2.75, height = 2, family = "serif")

    par(mar = c(0, 0, 0, 0))
    plot.new()

    text(0.175, 0.9, "Hyper", cex = 0.9, font = 2, srt = 25)
    text(0.35, 0.9, "Hypo", cex = 0.9, font = 2, srt = 25)

    first_4_TE <- function(x) {
        legend(x.cord - x, y.cord,
            legend = rep("", 6),
            fill = c("white", "white", "white", "white", "white", "#afad43"),
            border = c("white", "white", "white", "white", "white"),
            density = c(NA, NA, NA, NA, NA, 30), angle = c(NA, NA, NA, NA, NA, 30),
            bty = "n", cex = 1.1
        )
    }
    first_4_TE(0.045)
    first_4_TE(0.09)
    first_4_TE(0.135)


    legend(x.cord, y.cord,
        legend = c("CG", "CHG", "CHH", "", "CDS", "TEs"),
        fill = c("#009E73", "#0072B2", "#CC79A7", "white", "#000000", "#afad43"),
        border = c("#009E73", "#0072B2", "#CC79A7", "white", "#000000", "white"),
        density = c(30, 30, 30, NA, NA, 30), angle = c(30, 30, 30, NA, NA, 30),
        bty = "n", cex = 1.1, x.intersp = 1.5
    )

    legend(x.cord - 0.175, y.cord,
        legend = rep("", 6),
        fill = c("#009E73", "#0072B2", "#CC79A7", "white", "#000000", "#afad43"),
        border = c("#009E73", "#0072B2", "#CC79A7", "white", "#000000", "white"),
        density = c(NA, NA, NA, NA, NA, 30), angle = c(NA, NA, NA, NA, NA, 30),
        bty = "n", cex = 1.1
    )

    # line for CDS
    lines(c(0.15, 0.34), c(0.3, 0.3), lty = 1, lwd = 1.5, col = "#000000")

    # add box to 'TEs'
    rect(0.14, 0.1625, 0.3675, 0.215, border = "#afad43", lwd = 1.5)

    dev.off()
}

################################################################################

# build methylRawList from GRanges for each context
.jointGR_to_methylRawList <- function(gr, ids, context = "CpG", assembly = "custom") {
    keep <- switch(context,
        "CpG" = mcols(gr)$context %in% c("CG", "CpG"),
        "CHG" = mcols(gr)$context == "CHG",
        "CHH" = mcols(gr)$context == "CHH",
        stop("context must be CpG/CHG/CHH")
    )
    gr <- gr[keep]

    cn <- colnames(mcols(gr))
    mcol <- grep("^readsM[0-9]+$", cn, value = TRUE)
    idx <- as.integer(sub("^readsM", "", mcol))
    ncol <- paste0("readsN", idx)
    stopifnot(length(idx) == length(ids), all(ncol %in% cn))

    raws <- Map(function(mc, nc, sid) {
        m <- as.integer(mcols(gr)[[mc]])
        n <- as.integer(mcols(gr)[[nc]])
        cov <- m + n
        keep2 <- cov > 0
        df <- data.frame(
            chr      = as.character(seqnames(gr))[keep2],
            start    = start(gr)[keep2],
            end      = start(gr)[keep2],
            strand   = as.character(strand(gr))[keep2],
            coverage = cov[keep2],
            numCs    = m[keep2],
            numTs    = n[keep2]
        )
        new("methylRaw", df,
            sample.id = sid, assembly = assembly,
            context = context, resolution = "base"
        )
    }, mcol, ncol, ids)

    methylRawList(raws, treatment = rep(NA_integer_, length(raws)))
}

# returns a GRanges for DMPs_calc
methylkit_dmps_joint_as_gr <- function(gr_joint, condition, cntx = "CG",
                                       min_cov = 4,
                                       q_main = 0.05,
                                       diff_main = NULL,
                                       assembly = "custom") {
    stopifnot(length(unique(condition)) == 2)
    condition <- as.character(condition)

    # thresholds by context
    if (is.null(diff_main)) diff_main <- if (cntx == "CG") 20 else 10 # 20 else 10

    # IDs
    within_grp_idx <- ave(seq_along(condition), condition, FUN = seq_along)
    ids <- paste0(condition, within_grp_idx)

    # methylKit objects for the context
    ctx_mk <- if (cntx == "CG") "CpG" else cntx
    mr <- .jointGR_to_methylRawList(gr_joint, ids = ids, context = ctx_mk, assembly = assembly)

    grp <- factor(condition, levels = unique(condition))
    treat <- as.integer(grp == levels(grp)[2])
    mr@treatment <- treat

    # QC + unite
    mr <- filterByCoverage(mr, lo.count = min_cov, hi.perc = 99.9)
    mr <- normalizeCoverage(mr)

    min_per_group <- if (all(table(condition) >= 2)) 2L else 1L
    destrand <- (cntx == "CG")
    meth <- unite(mr, destrand = destrand, min.per.group = min_per_group)

    # stats (replicates -> overdispersion; else Fisher)
    diff_obj <- if (all(table(condition) >= 2)) {
        calculateDiffMeth(meth, overdispersion = "MN")
    } else {
        calculateDiffMeth(meth, test = "fast.fisher", overdispersion = "none")
    }

    # extract with main threshold; if empty, relax
    dmps_df <- getData(getMethylDiff(diff_obj, difference = diff_main, qvalue = q_main))
    if (!nrow(dmps_df)) {
        # Return empty GRanges with proper structure
        empty_gr <- GRanges()
        empty_counts_df <- data.frame(matrix(ncol = length(ids) * 2, nrow = 0))
        colnames(empty_counts_df) <- as.vector(rbind(
            paste0("readsM", seq_len(length(ids))),
            paste0("readsN", seq_len(length(ids)))
        ))

        mcols(empty_gr) <- DataFrame(
            context = factor(character(0), levels = cntx),
            empty_counts_df,
            trinucleotide_context = character(0),
            pValue = numeric(0),
            qValue = numeric(0),
            meth.diff = numeric(0),
            sumReadsM1 = integer(0),
            sumReadsN1 = integer(0),
            proportion1 = numeric(0),
            sumReadsM2 = integer(0),
            sumReadsN2 = integer(0),
            proportion2 = numeric(0),
            cytosinesCount = integer(0),
            direction = integer(0),
            regionType = character(0)
        )

        return(empty_gr)
    }

    # per-sample counts from methylBase (SAFE: use getData())
    meth_df <- getData(meth)
    mcn <- colnames(meth_df)
    numCs_cols <- grep("^numCs[0-9]+$", mcn, value = TRUE)
    numTs_cols <- grep("^numTs[0-9]+$", mcn, value = TRUE)
    s_idx <- order(as.integer(sub("^numCs", "", numCs_cols)))
    numCs_cols <- numCs_cols[s_idx]
    numTs_cols <- numTs_cols[s_idx]

    # merge by chr/start/end only (strand-agnostic)
    key_cols <- c("chr", "start", "end")
    merged <- merge(dmps_df,
        meth_df[, c(key_cols, numCs_cols, numTs_cols)],
        by = key_cols, all.x = TRUE, sort = FALSE
    )

    # keep 'merged' in the same row order as dmps_df
    key_dmps <- paste(dmps_df$chr, dmps_df$start, dmps_df$end)
    key_merged <- paste(merged$chr, merged$start, merged$end)
    merged <- merged[match(key_dmps, key_merged), , drop = FALSE]

    # rename per-sample columns to readsM#/readsN#
    readsM <- merged[, numCs_cols, drop = FALSE]
    readsN <- merged[, numTs_cols, drop = FALSE]
    colnames(readsM) <- paste0("readsM", seq_len(ncol(readsM)))
    colnames(readsN) <- paste0("readsN", seq_len(ncol(readsN)))
    counts_df <- cbind(readsM, readsN)[, as.vector(rbind(colnames(readsM), colnames(readsN))), drop = FALSE]

    # group sums / proportions
    g1 <- which(grp == levels(grp)[1])
    g2 <- which(grp == levels(grp)[2])
    sumReadsM1 <- rowSums(as.matrix(readsM[, g1, drop = FALSE]), na.rm = TRUE)
    sumReadsN1 <- rowSums(as.matrix(readsN[, g1, drop = FALSE]), na.rm = TRUE)
    sumReadsM2 <- rowSums(as.matrix(readsM[, g2, drop = FALSE]), na.rm = TRUE)
    sumReadsN2 <- rowSums(as.matrix(readsN[, g2, drop = FALSE]), na.rm = TRUE)
    proportion1 <- ifelse(sumReadsM1 + sumReadsN1 > 0, sumReadsM1 / (sumReadsM1 + sumReadsN1), NA_real_)
    proportion2 <- ifelse(sumReadsM2 + sumReadsN2 > 0, sumReadsM2 / (sumReadsM2 + sumReadsN2), NA_real_)

    # trinucleotide_context from joint GRanges (match by chr+start)
    gr_by_ctx <- if (cntx == "CG") {
        gr_joint[mcols(gr_joint)$context %in% c("CG", "CpG")]
    } else {
        gr_joint[mcols(gr_joint)$context == cntx]
    }
    key_ctx <- paste0(as.character(seqnames(gr_by_ctx)), ":", start(gr_by_ctx))
    key_out <- paste0(merged$chr, ":", merged$start)
    hit <- match(key_out, key_ctx)
    tri <- ifelse(is.na(hit), NA_character_,
        as.character(mcols(gr_by_ctx)$trinucleotide_context[hit])
    )

    # direction / regionType
    delta <- proportion2 - proportion1
    direction <- ifelse(delta > 0, 1L, ifelse(delta < 0, -1L, 0L))
    regionType <- ifelse(direction == 1L, "gain",
        ifelse(direction == -1L, "loss", "flat")
    )

    # final GRanges
    gr_out <- GRanges(
        seqnames = merged$chr,
        ranges = IRanges::IRanges(merged$start, merged$end),
        strand = "*"
    )

    mcols(gr_out) <- DataFrame(
        context = factor(cntx),
        counts_df,
        trinucleotide_context = tri,
        pValue = merged$pvalue,
        qValue = merged$qvalue, # add qvalue
        meth.diff = merged$meth.diff, # add methylation difference (group2 - group1)
        sumReadsM1 = as.integer(sumReadsM1),
        sumReadsN1 = as.integer(sumReadsN1),
        proportion1 = as.numeric(proportion1),
        sumReadsM2 = as.integer(sumReadsM2),
        sumReadsN2 = as.integer(sumReadsN2),
        proportion2 = as.numeric(proportion2),
        cytosinesCount = rep.int(1L, length(gr_out)),
        direction = as.integer(direction),
        regionType = as.character(regionType)
    )

    gr_out
}

################################################################################
