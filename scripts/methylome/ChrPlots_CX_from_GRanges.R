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
    main = ""
  )
  if (length(methylationProfiles) > 1) {
    for (i in 2:length(methylationProfiles)) {
      lines(pos, methylationProfiles[[i]]$Proportion,
        type = "o", col = col[i], lty = lty[i], pch = pch[i]
      )
    }
  }

  lines(pos, rep(y_mid, length(pos)), type = "o", col = "gray30", lty = 2, pch = 26)
  mtext(paste0("Chr ", chr.n), side = 1, line = 0, cex = 1)
}

###################################################################

ChrPlots_CX <- function(comparison_name, meth_var_list, meth_names, cntx, y_max = 1, y_mid = NULL, y_min = 0, italic_legend_names = TRUE, ylab_suffix = NULL, output_dir = ".") {
  ### color palette
  n.pal <- ifelse(length(meth_var_list) < 3, 3, length(meth_var_list))
  col_vec <- c("#00000090", "#bf682890", paste0(RColorBrewer::brewer.pal(n = n.pal, name = "Set1")[-5], "90"))

  ### ylab suffix - in addition to 'CNTX methylation'
  ### add 'ylab_suffix=(delta)' to get 'CG methylation (delta)'
  if (!is.null(ylab_suffix)) {
    ylab_suffix <- paste0(" ", ylab_suffix)
  }

  ### change seqnames for plot
  meth_vars_trimmed <- lapply(meth_var_list, function(m) renameSeqlevels(m, gsub("Chr", "", seqlevels(m))))

  ### change value column name
  meth_vars_trimmed <- lapply(meth_vars_trimmed, function(x) {
    names(mcols(x)) <- "Proportion"
    x
  })

  ### Low resolution profiles plot ###
  cat(paste0("\n", cntx, " ChrPlot..."))
  svg(paste0(output_dir, "/ChrPlot_", cntx, "_", comparison_name, ".svg"), width = 8, height = 2, family = "serif")
  par(mar = c(1, 4, 2, 0))
  par(fig = c(0, 2, 0, 10) / 10)
  plot(runif(10), runif(10),
    xlim = c(0, 0.01),
    ylim = c(y_min, y_max), axes = FALSE, type = "n", ylab = paste0(cntx, " methylation", ylab_suffix), xlab = ""
  )
  y_mid <- ifelse(is.null(y_mid), (y_max + y_min) / 2, y_mid)
  axis(2, c(y_min, y_mid, y_max), lty = 1, labels = c(y_min, y_mid, y_max))
  par(new = T)

  i <- 1
  u <- (10 - i) / 6
  for (chr_number in 1:5) {
    par(mar = c(1, 0, 2, 0))
    par(fig = c(i, i + u, 0, 10) / 10)
    chromosome_plot(meth_vars_trimmed, meth_names, chr_number, y_max, y_mid, y_min, col_vec)
    par(new = T)
    i <- i + u
  }
  par(mar = c(1, 0, 2, 0))
  par(fig = c(8, 10, 0, 10) / 10)
  legend("topright",
    legend = meth_names,
    text.font = ifelse(italic_legend_names, 3, 1),
    lty = 0, col = col_vec,
    pch = 16, bty = "n"
  )
  cat(" done\n")
  dev.off()
}
