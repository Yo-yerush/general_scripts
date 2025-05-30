chromosome_plot <- function(profile_vars, profile_names, chr.n, y_max, y_mid, y_min, col_vec, cntx) {
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
    if (cntx != "TE") {
        axis(1, at = c(min(pos), max(pos)), labels = FALSE, col = "gray40", tck = 0.025)
    } else {
        # color bellow the lines
        polygon(c(pos, rev(pos)),
        c(methylationProfiles[[1]]$Proportion, rep(y_min, length(pos))),
        col = adjustcolor(col[1], alpha.f = 0.2), border = NA)
        axis(1, at = c(min(pos), max(pos)), labels = FALSE, col = "gray20", tck = 0.1)
    }
    if (length(methylationProfiles) > 1) {
        for (i in 2:length(methylationProfiles)) {
            lines(pos, methylationProfiles[[i]]$Proportion,
                type = "o", col = col[i], lty = lty[i], pch = pch[i]
            )
        }
    }

    if (cntx != "TE") {
        lines(pos, rep(0, length(pos)), type = "o", col = "gray30", lty = 2, pch = 26)
    } # else {
    #   lines(pos, rep(0, length(pos)), type = "o", col = "black", lty = 1, pch = 26)
    # }

    # if (cntx == "CHH") {
    #    mtext(paste0("Chr ", chr.n), side = 1, line = 0, cex = 1)
    # }
}

###################################################################

ChrPlots_CX_all <- function(
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
    y_title_cex = 1,
    TE_as_gr = "tair10") {
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

    if (is.null(TE_as_gr)) {
        lay <- cbind(
            matrix(1:24, nrow = 4, byrow = TRUE),
            rep(25, 4) # legend in column 7, spans all rows
        )

        layout(lay,
            widths  = c(0.375, chr_len, 0.5),
            heights = c(1, 1, 1, 0.25) # last row for Chrs
        )
    } else {
        lay <- cbind(
            matrix(1:30, nrow = 5, byrow = TRUE),
            rep(31, 5)
        )

        layout(lay,
            widths  = c(0.375, chr_len, 0.5),
            heights = c(1, 1, 1, 0.35, 0.30)
        )
    }


    for (cntx in c("CG", "CHG", "CHH", "TE")) {
        meth_vars_context <- lapply(meth_vars_trimmed, function(inner_list) inner_list[[tolower(cntx)]])

        if (cntx != "TE") {
            ## y-axis
            y_title <- paste0(cntx, " methylation", ylab_suffix)
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
                adj = ifelse(y_min_cntx == 0, 0.5, 0), #
                line = 0.65,
                col = "gray25",
                cex = 0.55
            )
            mtext(
                side = 2,
                text = y_mid_cntx,
                at = y_mid_cntx,
                adj = 0.5, #
                line = 0.65,
                col = "gray25",
                cex = 0.55
            )
            mtext(
                side = 2,
                text = y_max_cntx,
                at = y_max_cntx,
                adj = ifelse(y_max_cntx == 0, 0.5, 0.9), #
                line = 0.65,
                col = "gray25",
                cex = 0.55
            )

            ## chromosome
            par(mar = c(0, 0, 0, 0)) # c(1, 0, 2, 0))
            for (chr in 1:5) {
                chromosome_plot(meth_vars_context, meth_names, chr, y_max_cntx, y_mid_cntx, y_min_cntx, col_vec, cntx)
            }
        } else {
            if (!is.null(TE_as_gr)) {
                if (TE_as_gr == "tair10") {
                    te_plot_conf(upload_te())
                } else {
                    te_plot_conf(TE_as_gr)
                }
            }
        }
    }

    ## x-axis - chromosome
    par(mar = c(1, 0, 0, 0))
    plot.new()
    for (chr in 1:5) {
        plot.new()
        mtext(
            side = 1,
            text = paste0("Chr ", chr),
            line = -0.5,
            at = 0.5,
            adj = 0.5,
            col = "gray25"
        )
    }

    ## legend
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("top",
        legend = meth_names,
        text.font = ifelse(italic_legend_names, 3, 1),
        col = sub(".{2}$", "", col_vec), # remove the last two characters (e.g., "90")
        lty = 1, bty = "n", cex = 1.2, lwd = 2
    )
}

###################################################################

upload_te <- function() {
    te_file <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_160525/annotation_files/TAIR10_Transposable_Elements.txt"
    source("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_160525/scripts/trimm_and_rename_seq.R")
    source("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/methionine/Methylome.At_paper/files_160525/scripts/edit_TE_file.R")

    TE_4_dens <- edit_TE_file(read.csv(te_file, sep = "\t")) %>%
        circlize::genomicDensity(window.size = 0.5e6) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    return(TE_4_dens)
}

te_plot_conf <- function(te_gr) {
    te_gr <- renameSeqlevels(te_gr, gsub("Chr|chr|chromosome", "", seqlevels(te_gr)))
    names(mcols(te_gr)) <- "Proportion"

    y_max_te <- 1 # max(te_gr$Proportion)
    y_mid_te <- 0.5 # y_max_te / 2
    y_min_te <- 0

    par(mar = c(0, 4, 0, 0))
    plot(NA, NA,
        xlim = c(0, 1), ylim = c(y_min_te, y_max_te),
        axes = FALSE, xlab = "",
        ylab = "TEs", cex.lab = 1.2
    )

    axis(2, at = c(y_min_te, y_max_te), labels = FALSE)
    mtext(
        side = 2, text = c("", ""), # c(y_min_te, y_mid_te, round(y_max_te, 2)),
        at = c(y_min_te, y_max_te) # ,
        # line = 0.65, col = "gray25", cex = 0.7, adj = c(0, .5, .9)
    )

    par(mar = c(0, 0, 0, 0))
    te_list <- list(te_gr)
    te_names <- "TE"

    for (chr in 1:5) {
        chromosome_plot(te_list, te_names,
            chr, y_max_te, y_mid_te, y_min_te,
            col_vec = "#55555590", cntx = "TE"
        )
    }
}
