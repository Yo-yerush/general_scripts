####### trimm 'ChrM' and 'ChrC'
trimm_Chr <- function(gr_obj) {
  remove_seqnames = c("NC_000932.1","NC_037304.1")
  gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
  seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
  return(sort(gr_obj))
}

####### 'DMRcaller' functions
.sumReadsM <- function(methylationData){
  return(sum(methylationData$readsM))
}

.sumReadsN <- function(methylationData){
  return(sum(methylationData$readsN))
}

.analyseReadsInsideRegionsOneSample <- function(methylationData, regions){
  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- S4Vectors::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))
  
  regions$sumReadsM <- rep(0, times=length(regions))
  regions$sumReadsN <- rep(0, times=length(regions))    
  regions$Proportion <- rep(0, times=length(regions))        
  
  
  regions$sumReadsM[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsM)
  regions$sumReadsN[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsN)  
  
  regions$Proportion[regionsIndexes] <- regions$sumReadsM[regionsIndexes]/regions$sumReadsN[regionsIndexes]      
  return(regions)
}


####### yo modified 'DMRcaller' scripts
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
  } else{
    contextMethylationData <- localMethylationData[localMethylationData$trinucleotide_context %in% 
                                                     context]
  }
  
  rm(localMethylationData)
  
  seqs = seq(minPos, maxPos - windowSize, windowSize)
  ranges <- GRanges(seqname, IRanges(seqs, seqs + windowSize - 
                                       1))
  
  
  ranges <- .analyseReadsInsideRegionsOneSample(contextMethylationData, 
                                                ranges)
  ranges$context <- paste(context, collapse = "_")
  return(ranges)
}

chromosome_plot <- function(delta_profile, #var1_profile, var2_profile,
                            chr.n, chr.name, cntx,
                            difference=F, y.new.scale=F, y_cntx_max=NULL, y_cntx_min=NULL) {
  
  # if not delta
  #var1_profile <- var1_profile[seqnames(var1_profile) %in% chr.name]
  #var2_profile <- var2_profile[seqnames(var2_profile) %in% chr.name]
  
  delta_profile <- delta_profile[seqnames(delta_profile) %in% chr.name]
  
  
  if (chr.n == 1) {ylab = "methylation"} else {ylab=''}
  
  col = brewer.pal(n = length(delta_profile), name = 'Set1')
  pch = rep(26,length(delta_profile))
  lty = rep(1,length(delta_profile))
  
  
  # if not delta
  #col_control = brewer.pal(n = length(delta_profile), name = 'Pastel1')
  #lty_control = rep(2,length(delta_profile))
  
  
  if (y.new.scale) {
    ymax = y_cntx_max
    ymin = y_cntx_min
  } else {
    ymax = 1
    ymin = 0
  }
  
  pos = (start(delta_profile[[1]]) + end(delta_profile[[1]]))/2
  # lines for treatment var
  plot(pos, delta_profile[[1]]$Proportion, type = "o", 
       ylim = c(ymin, ymax), xlab = "genomic coordinate", ylab = ylab, 
       col = col[1], pch = pch[1], lty = lty[1], 
       yaxt = "n", xaxt = "n", 
       main = ""
  )
  
  
  # lines for the rest of treatment vars
  for (i in 2:length(delta_profile)) {
    lines(pos, delta_profile[[i]]$Proportion, 
          type = "o", col = col[i], lty = lty[i], pch = pch[i])
  }
  
  # if not delta
  # dashed lines for control vars
  #for (i in 1:length(var1_profile)) {
  #  lines(pos, var1_profile[[i]]$Proportion, 
  #        type = "o", col = col_control[i], lty = lty_control[i], pch = pch[i])
  #}
  
  
  # line at 0 (yaxis)
  lines(pos, rep(0,length(delta_profile[[1]])), 
        type = "o", col = "gray30", lty = 2, pch = 26)
  
  mtext(paste0("Chr ",chr.n), side = 1, line = 0, cex = 1)#adj = 0,
  
}