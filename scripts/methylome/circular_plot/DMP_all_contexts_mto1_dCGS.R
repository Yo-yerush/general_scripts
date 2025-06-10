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

  
  mto1_samples = 9:11
  #dcgs_samples = c(20,22,23)
  #ev_samples = c(24:26)
  
  ############# CX report #############
  wt_cx = paste0("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S",18:19,"/methylation_extractor/S",18:19,"_R1_bismark_bt2_pe.CX_report.txt")
  mto1_cx = paste0("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S",mto1_samples,"/methylation_extractor/S",mto1_samples,"_R1_bismark_bt2_pe.CX_report.txt")
  
  dcgs_cx = c("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S20/methylation_extractor/V350126508_L04_ARAedutH012200-1_1_bismark_bt2_pe.CX_report.txt",
              "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S22/methylation_extractor/V350126508_L04_ARAedutH012201-2_1_bismark_bt2_pe.CX_report.txt",
              "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S23/methylation_extractor/V350126508_L04_ARAedutH012202-3_1_bismark_bt2_pe.CX_report.txt")
  ev_cx = c("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S24/methylation_extractor/V350126648_L01_ARAedutH012203-4_1_bismark_bt2_pe.CX_report.txt",
            "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S25/methylation_extractor/V350126648_L01_ARAedutH012204-13_1_bismark_bt2_pe.CX_report.txt",
            "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S26/methylation_extractor/V350126648_L01_ARAedutH012205-14_1_bismark_bt2_pe.CX_report.txt")
  
  wt_methylationData <- mclapply(wt_cx, readBismark, mc.cores = 2)
  mto1_methylationData <- mclapply(mto1_cx, readBismark, mc.cores = 3)
  dcgs_methylationData <- mclapply(dcgs_cx, readBismark, mc.cores = 3)
  ev_methylationData <- mclapply(ev_cx, readBismark, mc.cores = 3)
  
  wt_pool <- poolMethylationDatasets(GRangesList(wt_methylationData))
  mto1_pool <- poolMethylationDatasets(GRangesList(mto1_methylationData))
  dcgs_pool <- poolMethylationDatasets(GRangesList(dcgs_methylationData))
  ev_pool <- poolMethylationDatasets(GRangesList(ev_methylationData))
  
  wt_pool = trimm_Chr(wt_pool)
  mto1_pool = trimm_Chr(mto1_pool)
  dcgs_pool = trimm_Chr(dcgs_pool)
  ev_pool = trimm_Chr(ev_pool)
  
  
  delta_list_mto1 = list(CG=NA, CHG=NA, CHH=NA)
  delta_list_dcgs = list(CG=NA, CHG=NA, CHH=NA)
  
  for (context in c("CG","CHG","CHH")) {
    
    profile_wt = GRanges()
    profile_mto1 = GRanges()
    profile_dcgs = GRanges()
    profile_ev = GRanges()
    
    for (chr.n in as.character(wt_pool@seqnames@values)) {
      chr.vec = wt_pool[wt_pool@seqnames == chr.n]
      regions <- GRanges(seqnames = Rle(chr.n), ranges = IRanges(chr.vec@ranges@start[1],chr.vec@ranges@start[length(chr.vec@ranges@start)]))
      profile_wt <- c(profile_wt, computeMethylationProfile(wt_pool, regions, windowSize = 150000, context = context))
      profile_mto1 <- c(profile_mto1, computeMethylationProfile(mto1_pool, regions, windowSize = 150000, context = context))
      profile_dcgs <- c(profile_dcgs, computeMethylationProfile(dcgs_pool, regions, windowSize = 150000, context = context))
      profile_ev <- c(profile_ev, computeMethylationProfile(ev_pool, regions, windowSize = 150000, context = context))
    }
    methylationProfiles_mto1 = profile_wt[,0]
    methylationProfiles_dcgs = profile_wt[,0]
    
    methylationProfiles_mto1$Proportion = profile_mto1$Proportion - profile_wt$Proportion
    methylationProfiles_dcgs$Proportion = profile_dcgs$Proportion - profile_ev$Proportion
    
    delta_df_mto1 = as.data.frame(methylationProfiles_mto1)[,c("seqnames","start","end","Proportion")]
    delta_df_dcgs = as.data.frame(methylationProfiles_dcgs)[,c("seqnames","start","end","Proportion")]
    #delta_pos  = (delta_df$start + delta_df$end) / 2
    
    delta_df_dcgs[delta_df_dcgs < -0.025] = -0.025 # for the outlier one in Chr 3 around 14Mbp
    
    delta_list_mto1[[context]] = delta_df_mto1
    delta_list_dcgs[[context]] = delta_df_dcgs
  }
  
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

  #####################################
  #y_ranges = list(CG = range(c(delta_list_mto1[["CG"]]$Proportion, delta_list_dcgs[["CG"]]$Proportion), na.rm = TRUE),
  #                CHG = range(c(delta_list_mto1[["CHG"]]$Proportion, delta_list_dcgs[["CHG"]]$Proportion), na.rm = TRUE),
  #                CHH = range(c(delta_list_mto1[["CHH"]]$Proportion, delta_list_dcgs[["CHH"]]$Proportion), na.rm = TRUE))
  ############# the plot #############
  svg(paste0("/home/yoyerush/yo/methylome_pipeline/circular_plots/DMP_all_contexts/mto1_dcgs_all_contexts_circular_plots.svg"), width = 4.25, height = 4.25, family = "serif")
  circos.genomicInitialize(as.data.frame(gff3.trimmed)[,1:3], sector.names = paste0("Chr ",1:5), axis.labels.cex = 0.45, labels.cex = 1.35)
  
  i=2
  for (cntx_plot in c("CG", "CHG", "CHH")) {
    
    plot_df = cbind(delta_list_mto1[[cntx_plot]], delta_list_dcgs[[cntx_plot]][,4])
    names(plot_df)[4:5] = c("P_mto1", "P_dcgs")
    
    circos.genomicTrackPlotRegion(plot_df, bg.col = "#fafcff", bg.border = NA, panel.fun = function(region, value, ...) {
      ## draw 'y = 0' line
      xlim = get.cell.meta.data("xlim")
      circos.lines(xlim, c(0, 0), col = "gray", lwd = 1, lty = 2) # y = 0
      
      circos.genomicLines(region, value[["P_mto1"]], col = "black")
      circos.genomicLines(region, value[["P_dcgs"]], col = "#bf6828")
    }, track.height = 0.23, track.margin = c(0,0))

    i=i+1
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
