library(DMRcaller)
library(rtracklayer)
library(dplyr)
library(circlize)
library(parallel)


trimm_Chr <- function(gr_obj) {
  remove_seqnames = c("NC_000932.1","NC_037304.1")
  gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
  seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
  return(sort(gr_obj))
}

####### DMRcaller functions
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

####### annotation files
{
  ############# gff3 file #############
  gff3 = import.gff3("/home/yoyerush/yo/TAIR10.1/GTF_file/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
  gff3.trimmed = trimm_Chr(gff3)
  
  ############# TE file #############
  TE = read.csv("/home/yoyerush/yo/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt", sep = "\t")
  TE_4_dens = TE[,c(1,3,4)]
  TE_4_dens$Transposon_Name = gsub(paste0("AT1TE.*"), "NC_003070.9",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT2TE.*"), "NC_003071.7",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT3TE.*"), "NC_003074.8",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT4TE.*"), "NC_003075.7",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT5TE.*"), "NC_003076.8",TE_4_dens$Transposon_Name)
  
}

mto1_samples = 9:11

context_list = list(CG = c("CGA","CGT","CGC","CGG"),
                    CHG = c("CAG", "CTG", "CCG"),
                    CHH = c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC"))

############# CX report #############
wt_cx = paste0("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S",18:19,"/methylation_extractor/S",18:19,"_R1_bismark_bt2_pe.CX_report.txt")
mto1_cx = paste0("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S",mto1_samples,"/methylation_extractor/S",mto1_samples,"_R1_bismark_bt2_pe.CX_report.txt")

wt_methylationData <- mclapply(wt_cx, readBismark, mc.cores = 2)
mto1_methylationData <- mclapply(mto1_cx, readBismark, mc.cores = 3)

wt_pool <- poolMethylationDatasets(GRangesList(wt_methylationData))
mto1_pool <- poolMethylationDatasets(GRangesList(mto1_methylationData))

wt_pool = trimm_Chr(wt_pool)
mto1_pool = trimm_Chr(mto1_pool)

##################################################

context_proportion_list_wt = list()
context_proportion_list_mto1 = list()

for (context_name in names(context_list)) {

  
  subContext_list_wt = list()
  subContext_list_mto1 = list()
  
  for (sub_context in context_list[[context_name]]) {
    
    profile_wt = GRanges()
    profile_mto1 = GRanges()
    
    for (chr.n in as.character(wt_pool@seqnames@values)) {
      chr.vec = wt_pool[wt_pool@seqnames == chr.n]
      regions <- GRanges(seqnames = Rle(chr.n), ranges = IRanges(chr.vec@ranges@start[1],chr.vec@ranges@start[length(chr.vec@ranges@start)]))
      
      profile_wt = c(profile_wt, MethylationProfile(wt_pool, regions, windowSize = 150000, context = sub_context))
      profile_mto1 = c(profile_mto1, MethylationProfile(mto1_pool, regions, windowSize = 150000, context = sub_context))
      
    }
    
    subContext_list_wt[[sub_context]] = as.data.frame(profile_wt)[,c("seqnames","start","end","Proportion")]
    subContext_list_mto1[[sub_context]] = as.data.frame(profile_mto1)[,c("seqnames","start","end","Proportion")]
    
    ###### delta
    #methylationProfiles_mto1 = profile_wt[,0]
    #methylationProfiles_mto1$Proportion = profile_mto1$Proportion - profile_wt$Proportion
    #delta_df_mto1 = as.data.frame(methylationProfiles_mto1)[,c("seqnames","start","end","Proportion")]
    #delta_pos  = (delta_df$start + delta_df$end) / 2
    #delta_list_mto1[[context]] = delta_df_mto1
  }

  {
  ################################################################
  # create data-frame with all sub-context (in columns)
  create_df_all_sunCNTX <- function(subContext_list) {
    proportion_list = subContext_list[[1]]
    iii = 2
    while (iii <= length(subContext_list)) {
      proportion_list = cbind(proportion_list,
                              subContext_list[[iii]][,4])
      iii=iii+1
    }
    return(proportion_list)
  }
  
  #context_proportion_list_wt[[context_name]] = create_df_all_sunCNTX(subContext_list_wt)
  #context_proportion_list_mto1[[context_name]] = create_df_all_sunCNTX(subContext_list_mto1)
  ################################################################
}
#}

############# the plot #############
svg(paste0("/home/yoyerush/yo/methylome_pipeline/circular_plots/DMP_all_contexts/mto1_sub_contexts_circular_plots.svg"), width = 4.25, height = 4.25, family = "serif")
circos.genomicInitialize(as.data.frame(gff3.trimmed)[,1:3], sector.names = paste0("Chr ",1:5), axis.labels.cex = 0.45, labels.cex = 1.35)


#for (cntx_plot in c("CG", "CHG", "CHH")) {
  
  plot_df = cbind(delta_list_mto1[[cntx_plot]], delta_list_mto3[[cntx_plot]][,4])
  names(plot_df)[4:5] = c("P_mto1", "P_mto3")
  
  circos.genomicTrackPlotRegion(plot_df, bg.col = "#fafcff", bg.border = NA, panel.fun = function(region, value, ...) {
    ## draw 'y = 0' line
    #xlim = get.cell.meta.data("xlim")
    #circos.lines(xlim, c(0, 0), col = "gray", lwd = 1, lty = 2) # y = 0
    
    
    for (variable in vector) {
      circos.genomicLines(region, value[["P_mto1"]], col = "black")
      circos.genomicLines(region, value[["P_mto3"]], col = "#4b63d6")
    }

    
  }, track.height = 0.23, track.margin = c(0,0))
  
  
}

circos.genomicDensity(TE_4_dens[1:3],
                      bg.col = "#fafcff", bg.border = NA, count_by = "number",
                      col ="#fcba0320", border = T, track.height = 0.115, track.margin = c(0,0))

circos.clear()
dev.off()

#### legends ####
if (F) {
  counts_norm_hyper = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  counts_norm_hypo = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  rndm_norm_hyper = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0,10))
  rndm_norm_hypo = c(rep(0,10), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
  
  legend_titles_circular = c("CG","CHG","CHH")
  tracks_total_size = 40
  legend_tracks_pos = as.character((tracks_total_size/4)+1)
  text_tracks_pos = as.character((tracks_total_size/4))
  
  svg(paste0("/home/yoyerush/yo/methylome_pipeline/circular_plots/DMP_all_contexts/legends.svg"), width = 6, height = 4.5, family = "serif")
  #svg(paste0("P:/yonatan/methionine/circular_plot_res/corr_TE_exp/legends.svg"), width = 6, height = 5.75, family = "serif")
  
  layout(matrix(1:2,ncol=1,nrow=2),heights = c(1,0.25))
  
  par(mar = c(0, 0.1, 0, 0.1) )
  circos.par("track.height" = 0.7, "canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1), "gap.degree" = 0, "clock.wise" = FALSE)
  circos.initialize(factors = as.character(1:tracks_total_size), xlim = c(0, 1)) 
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.1) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/99), rnorm(100,0.5, 0.1), area = F, border = T)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[1], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.1) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/99), rnorm(100,0.5, 0.1), area = F, border = T)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[2], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.1) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/99), rnorm(100,0.5, 0.1), area = F, border = T)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[3], facing = "downward", adj = c(0, 0.5))
  
  circos.clear()
  
  dev.off()
  
}
