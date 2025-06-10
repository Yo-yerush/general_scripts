library(DMRcaller)
library(rtracklayer)
library(dplyr)
library(parallel)
library(RColorBrewer)

#####
##### make the x axis size as the Chr size!



################################################################
## write all function as scripts in a directory called 'CX_sub_contexts_ChrPlot' and change this script file name...

trimm_Chr <- function(gr_obj) {
  remove_seqnames <- c("NC_000932.1", "NC_037304.1")
  gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
  seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
  return(sort(gr_obj))
}

####### DMRcaller functions
.sumReadsM <- function(methylationData) {
  return(sum(methylationData$readsM))
}

.sumReadsN <- function(methylationData) {
  return(sum(methylationData$readsN))
}

.analyseReadsInsideRegionsOneSample <- function(methylationData, regions) {
  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- S4Vectors::splitAsList(methylationData[queryHits(overlaps)], subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))

  regions$sumReadsM <- rep(0, times = length(regions))
  regions$sumReadsN <- rep(0, times = length(regions))
  regions$Proportion <- rep(0, times = length(regions))


  regions$sumReadsM[regionsIndexes] <- sapply(methylationDataContextList, .sumReadsM)
  regions$sumReadsN[regionsIndexes] <- sapply(methylationDataContextList, .sumReadsN)

  regions$Proportion[regionsIndexes] <- regions$sumReadsM[regionsIndexes] / regions$sumReadsN[regionsIndexes]
  return(regions)
}

MethylationProfile <- function(methylationData, region, windowSize, context) {
  seqname <- seqnames(region)
  minPos <- start(region)
  maxPos <- end(region)
  hits <- findOverlaps(methylationData, region)
  localMethylationData <- methylationData[queryHits(hits)]
  rm(methylationData)

  if (context == "CG" | context == "CHG" | context == "CHH") {
    contextMethylationData <- localMethylationData[localMethylationData$context %in%
      context]
  } else {
    contextMethylationData <- localMethylationData[localMethylationData$trinucleotide_context %in%
      context]
  }

  rm(localMethylationData)

  seqs <- seq(minPos, maxPos - windowSize, windowSize)
  ranges <- GRanges(seqname, IRanges(seqs, seqs + windowSize -
    1))


  ranges <- .analyseReadsInsideRegionsOneSample(
    contextMethylationData,
    ranges
  )
  ranges$context <- paste(context, collapse = "_")
  return(ranges)
}

chromosome_plot <- function(delta_profile, # var1_profile, var2_profile,
                            chr.n, chr.name, cntx, trnt,
                            difference = F, y.new.scale = F, y_cntx_max = NULL, y_cntx_min = NULL) {
  # if not delta
  # var1_profile <- var1_profile[seqnames(var1_profile) %in% chr.name]
  # var2_profile <- var2_profile[seqnames(var2_profile) %in% chr.name]

  delta_profile <- delta_profile[seqnames(delta_profile) %in% chr.name]


  if (chr.n == 1) {ylab <- "methylation"} else {ylab <- ""}

  col <- brewer.pal(n = length(delta_profile), name = "Set1")
  pch <- rep(26, length(delta_profile))
  lty <- rep(1, length(delta_profile))
  lwd = rep(2,length(delta_profile))

  # if not delta
  # col_control = brewer.pal(n = length(delta_profile), name = 'Pastel1')
  # lty_control = rep(2,length(delta_profile))


  if (y.new.scale) {
    ymax <- y_cntx_max
    ymin <- y_cntx_min
  } else {
    ymax <- 1
    ymin <- 0
  }

  pos <- (start(delta_profile[[1]]) + end(delta_profile[[1]])) / 2
  # lines for treatment var
  plot(pos, delta_profile[[1]]$Proportion,
    type = "o",
    ylim = c(ymin, ymax), xlab = "genomic coordinate", ylab = ylab,
    col = col[1], pch = pch[1], lty = lty[1],
    yaxt = "n", xaxt = "n",
    main = ""
  )


  # lines for the rest of treatment vars
  for (i in 2:length(delta_profile)) {
    lines(pos, delta_profile[[i]]$Proportion,
      type = "o", col = col[i], lty = lty[i], pch = pch[i]
    )
  }

  # if not delta
  # dashed lines for control vars
  # for (i in 1:length(var1_profile)) {
  #  lines(pos, var1_profile[[i]]$Proportion,
  #        type = "o", col = col_control[i], lty = lty_control[i], pch = pch[i])
  # }


  # line at 0 (yaxis)
  lines(pos, rep(0, length(delta_profile[[1]])),
    type = "o", col = "gray30", lty = 2, pch = 26
  )

  mtext(paste0("Chr ", chr.n), side = 1, line = 0, cex = 1) # adj = 0,

  # Adding text at the top left corner
  if (chr.n == 1) {
    text(
      x = par("usr")[1] + 0.125 * diff(par("usr")[1:2]), # Move a small fraction of the plot width to the right,
      y = par("usr")[4],
      labels = trnt,
      pos = 1, # Position text to the right of the specified coordinates
      # offset = 1,  # Add some space between the text and the top left corner
      cex = 0.9,
      font = 2,
      adj = c(0, 1)
    )
  }
}

################################################################
####### annotation files
{
  ############# gff3 file #############
  gff3 <- import.gff3("/home/yoyerush/yo/TAIR10.1/GTF_file/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
  gff3.trimmed <- trimm_Chr(gff3)
  chr_names <- unique(as.character(seqnames(gff3.trimmed)))

  ############# TE file #############
  TE <- read.csv("/home/yoyerush/yo/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt", sep = "\t")
  TE_4_dens <- TE[, c(1, 3, 4)]
  TE_4_dens$Transposon_Name <- gsub(paste0("AT1TE.*"), "NC_003070.9", TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name <- gsub(paste0("AT2TE.*"), "NC_003071.7", TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name <- gsub(paste0("AT3TE.*"), "NC_003074.8", TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name <- gsub(paste0("AT4TE.*"), "NC_003075.7", TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name <- gsub(paste0("AT5TE.*"), "NC_003076.8", TE_4_dens$Transposon_Name)
}

################################################################

for (treatment in c("mto1", "mto3", "dCGS", "sseHigh", "sseLow")) {
  cat("\n*", treatment, "*\n\n")

  ############# CX report #############

  if (treatment == "sseHigh") {
    trnt_file_samples <- "sseHigh"
    treatment <- "SSE_high"
  } else if (treatment == "sseLow") {
    trnt_file_samples <- "sseLow"
    treatment <- "SSE_low"
  }

  sample_df <- read.table(paste0("/home/yoyerush/yo/methylome_pipeline/Methylome.At_240325/Methylome.At/samples_table/samples_table_", trnt_file_samples, ".txt"))

  # Use mclapply to read files in parallel
  read_CX <- mclapply(sample_df[, 2], readBismark, mc.cores = nrow(sample_df))
  gr_list <- GRangesList(read_CX)

  ctrl_pool <- poolMethylationDatasets(gr_list[which(sample_df[, 1] != treatment)])
  trnt_pool <- poolMethylationDatasets(gr_list[which(sample_df[, 1] == treatment)])

  ctrl_pool <- trimm_Chr(ctrl_pool)
  trnt_pool <- trimm_Chr(trnt_pool)

  cat("\n")
  ################################################################

  context_list <- list(
    CG = c("CGA", "CGT", "CGC", "CGG"),
    CHG = c("CAG", "CTG", "CCG"),
    CHH = c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC")
  )


  for (context_name in names(context_list)) {
    cat(context_name, "sub-contexts ... ")

    subContext_list_ctrl <- GRangesList()
    subContext_list_trnt <- GRangesList()
    subContext_list_delta <- GRangesList()

    for (sub_context in context_list[[context_name]]) {
      profile_ctrl <- GRanges()
      profile_trnt <- GRanges()

      for (chr.n in as.character(ctrl_pool@seqnames@values)) {
        chr.vec <- ctrl_pool[ctrl_pool@seqnames == chr.n]
        regions <- GRanges(seqnames = Rle(chr.n), ranges = IRanges(chr.vec@ranges@start[1], chr.vec@ranges@start[length(chr.vec@ranges@start)]))

        profile_ctrl <- c(profile_ctrl, MethylationProfile(ctrl_pool, regions, windowSize = 500000, context = sub_context))
        profile_trnt <- c(profile_trnt, MethylationProfile(trnt_pool, regions, windowSize = 500000, context = sub_context))
      }


      subContext_list_delta[[sub_context]] <- profile_ctrl
      subContext_list_delta[[sub_context]]$Proportion <- (profile_trnt$Proportion - profile_ctrl$Proportion)

      # if its not delta
      # subContext_list_ctrl[[sub_context]] = profile_ctrl
      # subContext_list_trnt[[sub_context]] = profile_trnt



      cat(".")
    }

    ################################################################
    ############# the plot #############
    # if (context_name == "CG") {
    #  y_cntx = 1
    # } else if (context_name == "CHG") {
    #  y_cntx = 0.5
    # } else if (context_name == "CHH") {
    #  y_cntx = 0.2
    # }

    max_proportions <- sapply(subContext_list_delta, function(gr) {
      max(gr$Proportion)
    })

    min_proportions <- sapply(subContext_list_delta, function(gr) {
      min(gr$Proportion)
    })

    # Get the maximum value across all GRanges objects
    y_cntx_max <- max(max_proportions) * 1.1
    y_cntx_min <- min(min_proportions) * 1.1

    dir.create("/home/yoyerush/yo/methylome_pipeline/Chr_plots/delta_sub_contexts", showWarnings = F)

    svg(paste0("/home/yoyerush/yo/methylome_pipeline/Chr_plots/delta_sub_contexts/", treatment, "_sub.", context_name, "_delta_ChrPlot.svg"),
      width = 8, height = 2, family = "serif"
    )


    par(mar = c(1, 4, 2, 0))
    par(fig = c(0, 2, 0, 10) / 10)
    plot(runif(10), runif(10),
      xlim = c(0, 0.01),
      ylim = c(y_cntx_min, y_cntx_max), axes = FALSE, type = "n", ylab = paste(context_name, " methylation (Δ)"), xlab = ""
    )
    axis(2, c(y_cntx_min, y_cntx_max),
      lty = 1,
      labels = c(round(y_cntx_min, 4), round(y_cntx_max, 4))
    )
    par(new = T)

    i <- 1
    u <- (10 - i) / 6
    for (chr_number in 1:length(chr_names)) {
      par(mar = c(1, 0, 2, 0))
      par(fig = c(i, i + u, 0, 10) / 10)
      chromosome_plot(subContext_list_delta, # subContext_list_ctrl, subContext_list_trnt,
        chr_number, chr_names[chr_number], context_name, treatment,
        difference = F, y.new.scale = T, y_cntx_max = y_cntx_max, y_cntx_min = y_cntx_min
      )
      par(new = T)
      i <- i + u
    }

    ### legend
    par(mar = c(1, 0, 1, 0))
    par(fig = c(8, 10, 0, 10) / 10)
    legend("top",
      legend = context_list[[context_name]],
      lty = 0, col = brewer.pal(n = length(context_list[[context_name]]), name = "Set1"),
      pch = 15, bty = "n", cex = 0.75
    )

    dev.off()
    cat(" done\n")
  }
}
