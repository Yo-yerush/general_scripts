chromosome_plot <- function(profile_vars, profile_names, chr.n, y_max, y_mid, y_min, col_vec) {
    ### filter by chromosome
    profile_vars <- lapply(profile_vars, function(x) {
        x[x@seqnames == chr.n]
    })

    ### if the seqlevels are not the same, align them
    ### by taking the max length of each seqlevel
    maxLens <- Reduce(
        function(a, b) pmax(a, b, na.rm = TRUE),
        lapply(profile_vars, seqlengths)
    )
    common_si <- Seqinfo(names(maxLens), maxLens)
    profile_vars <- lapply(profile_vars, function(gr) {
        seqlevels(gr) <- seqlevels(common_si) # ensure same ordering
        seqinfo(gr) <- common_si # overwrite lengths
        gr
    })

    ### make it as gr list
    methylationProfiles <- GRangesList()
    for (i in seq_along(profile_vars)) {
        methylationProfiles[[profile_names[i]]] <- profile_vars[[i]]
    }

    col <- col_vec
    pch <- rep(26, length(methylationProfiles))
    lty <- rep(1, length(methylationProfiles))

    pos <- (start(methylationProfiles[[1]]) + end(methylationProfiles[[1]])) / 2
    plot(pos, methylationProfiles[[1]]$Proportion,
        type = "o",
        ylim = c(y_min, y_max), xlab = "", ylab = "",
        col = col[1], pch = pch[1], lty = lty[1],
        yaxt = "n", xaxt = "n",
        main = "", frame.plot = FALSE
    )
    axis(1, at = c(min(pos), max(pos)), labels = FALSE, col = "gray", tck = 0)

    if (length(methylationProfiles) > 1) {
        for (i in 2:length(methylationProfiles)) {
            lines(pos, methylationProfiles[[i]]$Proportion,
                type = "o", col = col[i], lty = lty[i], pch = pch[i]
            )
        }
    }

    lines(pos, rep(0, length(pos)), type = "o", col = "gray30", lty = 2, pch = 26)

    # if (cntx == "CHH") {
    #    mtext(paste0("Chr ", chr.n), side = 1, line = 0, cex = 1)
    # }
}

###################################################################

ChrPlots_CX_all <- function(
    comparison_name,
    meth_var_list,
    meth_names,
    y_max_cg = 1,
    y_max_chg = 0.5,
    y_max_chh = 0.2,
    y_mid_cg = NULL,
    y_mid_chg = NULL,
    y_mid_chh = NULL,
    y_min_cg = -1,
    y_min_chg = -0.5,
    y_min_chh = -0.2,
    italic_legend_names = TRUE,
    ylab_suffix = NULL,
    y_title_cex = 1.2,
    output_dir = ".") {
    ### color palette
    n.pal <- ifelse(length(meth_var_list) < 3, 3, length(meth_var_list))
    col_vec <- c("#00000090", "#bf682890", paste0(RColorBrewer::brewer.pal(n = n.pal, "Set1")[-5], "90"))

    ### y-mid edit
    y_mid_cg <- ifelse(is.null(y_mid_cg), (y_max_cg + y_min_cg) / 2, y_mid_cg)
    y_mid_chg <- ifelse(is.null(y_mid_chg), (y_max_chg + y_min_chg) / 2, y_mid_chg)
    y_mid_chh <- ifelse(is.null(y_mid_chh), (y_max_chh + y_min_chh) / 2, y_mid_chh)

    ### ylab suffix - in addition to 'CNTX methylation'
    ### add 'ylab_suffix=(delta)' to get 'CG methylation (delta)'
    if (!is.null(ylab_suffix)) {
        ylab_suffix <- paste0(" ", ylab_suffix)
    }

    ### change seqnames and value column name for plot
    meth_vars_trimmed <- lapply(meth_var_list, function(inner_list) {
        lapply(inner_list, function(m) {
            m <- renameSeqlevels(m, gsub("Chr", "", seqlevels(m)))
            names(mcols(m)) <- "Proportion"
            m
        })
    })

    ### normelize cheomosome panel to its size
    seq_len_vec <- as.numeric(seqlengths(meth_vars_trimmed[[1]][["chh"]]))
    chr_len <- seq_len_vec / max(seq_len_vec)

    ### Low resolution profiles plot ###
    ## column 1: y-axis strip
    ## columns 2-6: chromosomes 1-5
    ## column 7: legend
    svg(paste0(output_dir, "ChrPlot_", comparison_name, ".svg"), width = 7, height = 4, family = "serif")

    lay <- cbind(
        matrix(1:18, nrow = 3, byrow = TRUE),
        c(19, 19, 19) # legend
    )
    layout(lay,
        widths  = c(0.375, chr_len * 1, 0.5), # rep(1, 5)
        heights = rep(1, 4)
    )

    for (cntx in c("CG", "CHG", "CHH")) {
        meth_vars_context <- lapply(meth_vars_trimmed, function(inner_list) inner_list[[tolower(cntx)]])

        ## y-axis
        y_title <- paste0(cntx, " methylation ", ylab_suffix)
        y_max_cntx <- ifelse(cntx == "CG", y_max_cg, ifelse(cntx == "CHG", y_max_chg, y_max_chh))
        y_mid_cntx <- ifelse(cntx == "CG", y_mid_cg, ifelse(cntx == "CHG", y_mid_chg, y_mid_chh))
        y_min_cntx <- ifelse(cntx == "CG", y_min_cg, ifelse(cntx == "CHG", y_min_chg, y_min_chh))

        par(mar = c(0, 4, 0, 0)) # c(1, 4, 2, 0))
        plot(NA, NA,
            xlim = c(0, 1), ylim = c(y_min_cntx, y_max_cntx),
            axes = FALSE, xlab = "",
            ylab = y_title,
            cex.lab = y_title_cex
        )
        axis(2,
            at = c(y_min_cntx, y_mid_cntx, y_max_cntx),
            labels = FALSE # ,  c(paste0("       ", y_min_cntx), y_mid_cntx, paste0(y_max_cntx, "      ")),
            # col = "gray35", cex.axis = 1
        )
        mtext(
            side = 2,
            text = y_min_cntx,
            at = y_min_cntx,
            adj = 0, #
            line = 0.65,
            col = "gray25",
            cex = 0.7
        )
        mtext(
            side = 2,
            text = y_mid_cntx,
            at = y_mid_cntx,
            adj = 0.5, #
            line = 0.65,
            col = "gray25",
            cex = 0.7
        )
        mtext(
            side = 2,
            text = y_max_cntx,
            at = y_max_cntx,
            adj = 0.9, #
            line = 0.65,
            col = "gray25",
            cex = 0.7
        )

        ## chromosome
        par(mar = c(0, 0, 0, 0)) # c(1, 0, 2, 0))
        for (chr in 1:5) {
            chromosome_plot(meth_vars_context, meth_names, chr, y_max_cntx, y_mid_cntx, y_min_cntx, col_vec)
        }

    }

    ## legend
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("top",
        legend = meth_names,
        text.font = ifelse(italic_legend_names, 3, 1),
        lty = 1, col = col_vec, bty = "n", cex = 1.2, lwd = 3
    )

    dev.off()
}

###################################################################

upload_te <- function(te_file="C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_160525/annotation_files/TAIR10_Transposable_Elements.txt") {
    source("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_160525/scripts/trimm_and_rename_seq.R")
    source("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_160525/scripts/edit_TE_file.R")

    TE_4_dens <- edit_TE_file(read.csv(te_file, sep = "\t")) %>%
        circlize::genomicDensity(window.size = 1e8) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    TE_4_dens <- renameSeqlevels(TE_4_dens, gsub("Chr", "", seqlevels(TE_4_dens)))
    names(mcols(TE_4_dens)) <- "Proportion"

    return(list(TE_4_dens))
}
