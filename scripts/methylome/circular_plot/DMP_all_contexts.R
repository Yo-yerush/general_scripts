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

for (treatment in c("mto1","mto3")) {
  
  if (treatment == "mto1") {
    mto_samples = 9:11
  } else if (treatment == "mto3") {
    mto_samples = 12:14
  }
  ############# CX report #############
  wt_cx = paste0("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S",18:19,"/methylation_extractor/S",18:19,"_R1_bismark_bt2_pe.CX_report.txt")
  mto_cx = paste0("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S",mto_samples,"/methylation_extractor/S",mto_samples,"_R1_bismark_bt2_pe.CX_report.txt")
  wt_methylationData <- mclapply(wt_cx, readBismark, mc.cores = 2)
  mto_methylationData <- mclapply(mto_cx, readBismark, mc.cores = 3)
  wt_pool <- poolMethylationDatasets(GRangesList(wt_methylationData))
  mto_pool <- poolMethylationDatasets(GRangesList(mto_methylationData))
  wt_pool = trimm_Chr(wt_pool)
  mto_pool = trimm_Chr(mto_pool)
  
  
  delta_list = list(CG=NA, CHG=NA, CHH=NA)
  for (context in c("CG","CHG","CHH")) {
    
    profile_wt = GRanges()
    profile_mto = GRanges()
    for (chr.n in as.character(wt_pool@seqnames@values)) {
      chr.vec = wt_pool[wt_pool@seqnames == chr.n]
      regions <- GRanges(seqnames = Rle(chr.n), ranges = IRanges(chr.vec@ranges@start[1],chr.vec@ranges@start[length(chr.vec@ranges@start)]))
      profile_wt <- c(profile_wt, computeMethylationProfile(wt_pool, regions, windowSize = 150000, context = context))
      profile_mto <- c(profile_mto, computeMethylationProfile(mto_pool, regions, windowSize = 150000, context = context))
    }
    methylationProfiles = profile_wt[,0]
    methylationProfiles$Proportion = profile_mto$Proportion - profile_wt$Proportion
    delta_df = as.data.frame(methylationProfiles)[,c("seqnames","start","end","Proportion")]
    #delta_pos  = (delta_df$start + delta_df$end) / 2
    
    delta_list[[context]] = delta_df
  }
    
    ############# gff3 file #############
    gff3 = import.gff3("/home/yoyerush/yo/TAIR10.1/GTF_file/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
    gff3.trimmed = trimm_Chr(gff3)
    
    
    
    #####################################
    ############# the plot #############
    svg(paste0("/home/yoyerush/yo/methylome_pipeline/circular_plots/DMP_all_contexts/",treatment,"_all_contexts_circular_plots.svg"), width = 7.07, height = 7.40, family = "serif")
    circos.genomicInitialize(as.data.frame(gff3.trimmed)[,1:3], sector.names = paste0("Chr ",1:5), axis.labels.cex = 0.75, labels.cex = 1.5)
    
    circos.genomicTrackPlotRegion(delta_list[["CG"]], bg.col = "#fafcff",bg.border = NA, panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, col = "black")
    })
    
    circos.genomicTrackPlotRegion(delta_list[["CHG"]], bg.col = "#fafcff",bg.border = NA, panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, col = "black")
    })
    
    circos.genomicTrackPlotRegion(delta_list[["CHH"]], bg.col = "#fafcff",bg.border = NA, panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, col = "black")
    })
    
    circos.clear()
    dev.off()
}

#### legends ####
{
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
