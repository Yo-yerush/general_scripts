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
    delta_pos = (delta_df$start + delta_df$end) / 2
    #source("/home/yoyerush/yo/methylome_pipeline/Methylome.At/Methylome.At_scripts/ChrPlots_CX_fun.R")
    ############# gff3 file #############
    gff3 = import.gff3("/home/yoyerush/yo/TAIR10.1/GTF_file/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
    gff3.trimmed = trimm_Chr(gff3)
    
    transcripts = gff3.trimmed[which(gff3.trimmed$type == "transcript")]
    
    gff3_df = as.data.frame(gff3.trimmed)[,1:3]
    
    ############# TE file #############
    TE = read.csv("/home/yoyerush/yo/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt", sep = "\t")
    TE_4_dens = TE[,c(1,3,4)]
    TE_4_dens$Transposon_Name = gsub(paste0("AT1TE.*"), "NC_003070.9",TE_4_dens$Transposon_Name)
    TE_4_dens$Transposon_Name = gsub(paste0("AT2TE.*"), "NC_003071.7",TE_4_dens$Transposon_Name)
    TE_4_dens$Transposon_Name = gsub(paste0("AT3TE.*"), "NC_003074.8",TE_4_dens$Transposon_Name)
    TE_4_dens$Transposon_Name = gsub(paste0("AT4TE.*"), "NC_003075.7",TE_4_dens$Transposon_Name)
    TE_4_dens$Transposon_Name = gsub(paste0("AT5TE.*"), "NC_003076.8",TE_4_dens$Transposon_Name)
    
    ############# significant DMRs and RNAseq density #############
    RNA_file = read.csv(paste0("/home/yoyerush/yo/methylome_pipeline/circular_plots/RNAseq_res/all.transcripts.",treatment,"_vs_wt.DE.csv"))
    RNA_file = RNA_file[RNA_file$padj < 0.05, c("transcript_id","log2FoldChange")] %>% na.omit()
    names(RNA_file)[2] = "RNA_log2FC"
    
    transcript2merge = as.data.frame(transcripts)[,c("seqnames","start","end","ID")]
    names(transcript2merge)[4] = "transcript_id"
    transcript2merge$transcript_id = gsub("transcript:","",transcript2merge$transcript_id)
    RNA_transcriptsLoc = merge.data.frame(transcript2merge,RNA_file, by = "transcript_id")
    
    meth_file = read.csv(paste0("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/",treatment,"_vs_wt/DMRs_",context,"_",treatment,"_vs_wt.csv"))
    meth_file$log2FoldChange = log2(meth_file$proportion2 / meth_file$proportion1)
    meth_file = meth_file[,c("seqnames","start","end","log2FoldChange")]
    names(meth_file)[4] = "meth_log2FC"
    
    
    #####################################
    ############# the plot #############
    svg(paste0("/home/yoyerush/yo/methylome_pipeline/circular_plots/",treatment,"_",context,"_circular_plots.svg"), width = 7.07, height = 7.40, family = "serif")
    circos.genomicInitialize(as.data.frame(gff3.trimmed)[,1:3], sector.names = paste0("Chr ",1:5), axis.labels.cex = 0.75, labels.cex = 1.5)
    
    #if (treatment == "mto1" & context == "CHG") {
    #  circos.genomicDensity(list(meth_file[meth_file$meth_log2FC > 0,1:3],
    #                             meth_file[meth_file$meth_log2FC < 0,1:3]),
    #                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
    #                        col = c("#d96c6c","#6c96d9"))
    #} else if (treatment == "mto3" & context == "CHG") {
    #  circos.genomicDensity(list(meth_file[meth_file$meth_log2FC < 0,1:3],
    #                             meth_file[meth_file$meth_log2FC > 0,1:3]),
    #                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
    #                        col = c("#6c96d9","#d96c6c"))
    #} else if (treatment == "mto3" & context == "CG") {
    #  circos.genomicDensity(list(meth_file[meth_file$meth_log2FC < 0,1:3],
    #                             meth_file[meth_file$meth_log2FC > 0,1:3]),
    #                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
    #                        col = c("#6c96d9","#d96c6c"))
    #} else {
    #  circos.genomicDensity(list(meth_file[meth_file$meth_log2FC > 0,1:3],
    #                             meth_file[meth_file$meth_log2FC < 0,1:3]),
    #                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
    #                        col = c("#d96c6c","#6c96d9"))
    #}
    
      circos.genomicDensity(list(meth_file[meth_file$meth_log2FC > 0, 1:3],
                                 meth_file[meth_file$meth_log2FC < 0, 1:3]),
                            bg.col = "#fafcff", bg.border = NA, count_by = "number",
                            col = c("#FF000080","#304ed180"), border = T)
    
    #circos.genomicDensity(meth_file[meth_file$meth_log2FC > 0,2:4],
    #                      bg.col = "#fafcff", bg.border = NA, count_by = "number",
    #                      col = "#d96c6c", track.height = 0.1)
    #circos.genomicDensity(meth_file[meth_file$meth_log2FC < 0,2:4],
    #                      bg.col = "#fafcff", bg.border = NA, count_by = "number",
    #                      col = "#6c96d9", track.height = 0.1)
    
    circos.genomicTrackPlotRegion(delta_df, bg.col = "#fafcff",bg.border = NA, panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, col = "black")
    })
    
    #circos.genomicDensity(meth_file[meth_file$meth_log2FC > 0,2:4], count_by = "number", col = "#d96c6c")
    #circos.genomicDensity(meth_file[meth_file$meth_log2FC < 0,2:4], count_by = "number", col = "#6c96d9")
    circos.genomicDensity(RNA_transcriptsLoc[2:4], bg.col = "#fafcfc", bg.border = NA, count_by = "number", track.height = 0.175)
    circos.genomicDensity(TE_4_dens[1:3], bg.col = "#fafcfc", bg.border = NA, count_by = "number")
    circos.clear()
    dev.off()
  }
}

#### legends ####
{
  counts_norm_hyper = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  counts_norm_hypo = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  rndm_norm_hyper = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0,10))
  rndm_norm_hypo = c(rep(0,10), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
  
  legend_titles_circular = c("DMRs","?? Methylation","sig. Genes","TEs")
  tracks_total_size = 16
  legend_tracks_pos = as.character((tracks_total_size/4)+1)
  text_tracks_pos = as.character((tracks_total_size/4))
  
  svg(paste0("/home/yoyerush/yo/methylome_pipeline/circular_plots/legends.svg"), width = 6, height = 4.5, family = "serif")
  #svg(paste0("P:/yonatan/methionine/circular_plot_res/corr_TE_exp/legends.svg"), width = 6, height = 5.75, family = "serif")
  
  layout(matrix(1:2,ncol=1,nrow=2),heights = c(1,0.25))
  
  par(mar = c(0, 0.1, 0, 0.1) )
  circos.par("track.height" = 0.7, "canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1), "gap.degree" = 0, "clock.wise" = FALSE)
  circos.initialize(factors = as.character(1:tracks_total_size), xlim = c(0, 1)) 
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, area = T, border = T, col = "#FF000080")
  circos.lines(seq(0,1,1/(length(rndm_norm_hypo)-1)), rndm_norm_hypo, area = T, border = T, col = "#304ed180")
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[1], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.1) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/99), rnorm(100,0.5, 0.1), area = F, border = T)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[2], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/14), rnorm(15,0.5, 0.125), area = T, border = F)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[3], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/14), rnorm(15,0.5, 0.125), area = T, border = F)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[4], facing = "downward", adj = c(0, 0.5))
  
  circos.clear()
  
  par(mar = c(0, 0.1, 0, 0.1) )
  plot.new()
  legend(0.1,1,
         legend=c("Hyper","Hypo","Shared"),
         fill = c("#FF000095","#304ed195","#5e0d3d99"),
         bty="n",
         title=expression(bold("DMRs direction")))
  
  #plot.new()
  #  legend(0.2,1.05,
  #         legend=c(1,0,-1),
  #         fill = colorRampPalette(c("red","white", "blue"))(3),
  #        bty="n",
  #         title=expression(bold("Values")),
  #         title.adj=1)
  dev.off()
  
}
