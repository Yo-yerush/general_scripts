library(dplyr)
library(RColorBrewer)
library(rtracklayer)
library(circlize)

### GFF and TE files
{
  trimm_Chr <- function(gr_obj) {
    remove_seqnames = c("NC_000932.1","NC_037304.1")
    gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
    seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
    return(sort(gr_obj))
  }
  
  ############# gff3 fNULL############# gff3 file #############
  gff3 = import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
  gff3.trimmed = trimm_Chr(gff3)
  #genes = gff3.trimmed[which(gff3.trimmed$type == "gene")]
  
  transcripts = gff3.trimmed[which(gff3.trimmed$type == "transcript")]
  transcript2merge = as.data.frame(transcripts)[,c("seqnames","start","end","ID","gene")]
  names(transcript2merge)[4] = "transcript_id"
  transcript2merge$transcript_id = gsub("transcript:","",transcript2merge$transcript_id)
  
  
  ############# TE file #############
  TE = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/TAIR10_Transposable_Elements.txt", sep = "\t")
  TE_4_dens = TE[,c("Transposon_Name","Transposon_min_Start","Transposon_max_End","Transposon_Super_Family")]
  TE_4_dens$Transposon_Name = gsub(paste0("AT1TE.*"), "NC_003070.9",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT2TE.*"), "NC_003071.7",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT3TE.*"), "NC_003074.8",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT4TE.*"), "NC_003075.7",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT5TE.*"), "NC_003076.8",TE_4_dens$Transposon_Name)
} 

for (TE_SF_plot in c("gypsy_copia_L1","not_gypsy","LTR/Copia","LTR/Gypsy","LINE/L1")) {
  treatment = "mto1"
  
  cntx_file <- function(context) {
    ############# DMRs file
    dmrs_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/",treatment,"_vs_wt/DMRs_",context,"_",treatment,"_vs_wt.csv"))
    dmrs_file$log2FoldChange = log2(dmrs_file$proportion2 / dmrs_file$proportion1)
    dmrs_file = dmrs_file[,c("seqnames","start","end","log2FoldChange")]
    names(dmrs_file)[4] = "meth_log2FC"
    
    return(dmrs_file)
  }
  CG_file = cntx_file("CG")
  CHG_file = cntx_file("CHG")
  CHH_file = cntx_file("CHH")
  
  # bind this contexts data frames and make it as GenomicRanges object. each position make it 1 bp
  CG_gr = makeGRangesFromDataFrame(CG_file)
  CHG_gr = makeGRangesFromDataFrame(CHG_file)
  CHH_gr = makeGRangesFromDataFrame(CHH_file)
  
  #CG_gr = resize(CG_gr, width = 1)
  #CHG_gr = resize(CHG_gr, width = 1)
  #CHH_gr = resize(CHH_gr, width = 1)
  
  all_cntx_gr = c(CG_gr, CHG_gr, CHH_gr)
  all_cntx_gr = all_cntx_gr[!duplicated(all_cntx_gr)]
  
  
  ####### TE and DMRs overlay #######
  TE_gr = TE_4_dens
  names(TE_gr)[1:3] = c("seqnames","start","end")
  TE_gr$strand = "*"
  
  if (TE_SF_plot == "gypsy_copia_L1") {
    TE_gr = TE_gr[which(TE_gr$Transposon_Super_Family %in% c("LTR/Copia","LTR/Gypsy","LINE/L1")),]
    TE_gr$col = "#00000090"
    TE_gr$col[TE_gr$Transposon_Super_Family == "LTR/Copia"] = "#17870e90"
    TE_gr$col[TE_gr$Transposon_Super_Family == "LINE/L1"] = "#ab640290"
    
  } else if (TE_SF_plot == "not_gypsy") {
    TE_gr = TE_gr[which(TE_gr$Transposon_Super_Family %in% c("LTR/Copia","LTR/Gypsy","LINE/L1")),]
    TE_gr$col[TE_gr$Transposon_Super_Family == "LTR/Copia"] = "#17870e90"
    TE_gr$col[TE_gr$Transposon_Super_Family == "LINE/L1"] = "#ab640290"
  } else {
    TE_gr = TE_gr[TE_gr$Transposon_Super_Family == TE_SF_plot,]
    TE_gr$col = "#000000"
  }
  
  TE_gr = makeGRangesFromDataFrame(TE_gr, keep.extra.columns = T)
  
  m = findOverlaps(TE_gr, all_cntx_gr)
  overlap_TEs = TE_gr[queryHits(m)]
  overlap_TEs = resize(overlap_TEs, width = 1)
  overlap_TEs = overlap_TEs[!duplicated(overlap_TEs)]
  overlap_TEs_df = as.data.frame(overlap_TEs)
  
  dir.create("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/dmrs_all_cntx_TEs_overlapping_superFamilies/", showWarnings = F)
  dir.create(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/dmrs_all_cntx_TEs_overlapping_superFamilies/",treatment), showWarnings = F)
  #####################################
  ############# the plot #############
  svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/dmrs_all_cntx_TEs_overlapping_superFamilies/",
             treatment,"/",
             gsub("/",".",TE_SF_plot),
             "_TEs_overlapping_superFamilies_DMRs_",treatment,".svg"),
      width = 3.25, height = 3.25, family = "serif")
  
  circos.genomicInitialize(as.data.frame(gff3.trimmed)[,1:3], sector.names = paste0("Chr ",1:5), axis.labels.cex = 0.4, labels.cex = 1.35)
  
  circos.genomicDensity(list(CG_file[CG_file$meth_log2FC > 0, 1:3],
                             CG_file[CG_file$meth_log2FC < 0, 1:3]),
                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
                        col = c("#FF000080","#304ed180"), border = T, track.height = 0.16, track.margin=c(0, 0))
  
  circos.genomicDensity(list(CHG_file[CHG_file$meth_log2FC > 0, 1:3],
                             CHG_file[CHG_file$meth_log2FC < 0, 1:3]),
                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
                        col = c("#FF000080","#304ed180"), border = T, track.height = 0.16, track.margin=c(0, 0))
  
  circos.genomicDensity(list(CHH_file[CHH_file$meth_log2FC > 0, 1:3],
                             CHH_file[CHH_file$meth_log2FC < 0, 1:3]),
                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
                        col = c("#FF000080","#304ed180"), border = T, track.height = 0.16, track.margin=c(0, 0))
  
  circos.genomicTrackPlotRegion(overlap_TEs_df[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, panel.fun = function(region, value, ...) {
    region_colors <- overlap_TEs_df$col[which(overlap_TEs_df[,1] == get.current.chromosome())]
    circos.genomicLines(region, 1, type = "h", col = region_colors, lwd = 0.4)
  }, track.height = 0.085, track.margin=c(0, 0))
  
  circos.genomicDensity(TE_4_dens[1:3],
                        bg.col = "#fafcff", bg.border = NA, count_by = "number",
                        col ="#fcba0330", border = T, track.height = 0.2, track.margin=c(0, 0))
  
  circos.clear()
  
  dev.off()
  
  cat(">")
}
############################################################################################################
############################################################################################################

# legends
if (F) {
  counts_norm_hyper = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  counts_norm_hypo = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  rndm_norm_hyper = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0,10))
  rndm_norm_hypo = c(rep(0,10), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
  
  #rndm_norm_TE = c(rep(0,10),(counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0,10))
  rndm_norm_TE = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
  #rndm_norm_TE = c(rep(0,10), rndm_norm_hypo[-c(1:10, (length(rndm_norm_hypo)-10):length(rndm_norm_hypo))], rep(0,10))
  
  rndm_norm_genes = (1-rndm_norm_TE)/1.4
  #rndm_norm_genes = 0.6-rndm_norm_TE
  #rndm_norm_genes[rndm_norm_genes<0.12] = 0.12
  
  rndm = runif(12, min = 0, max = 1)
  legend_titles_circular = c("DMRs","  Group sig. genes","Methylation values","Expression values","  Overlapping TEs","Genome annotation")
  tracks_total_size = 16
  legend_tracks_pos = as.character((tracks_total_size/4)+1)
  text_tracks_pos = as.character((tracks_total_size/4))
  
  colfunc <- colorRampPalette(c("red", "white", "blue"))
  colfunc_vec_1 = colfunc(16)[-c((round(16/2)-1):(16/2), ((16/2)+1):((16/2)+1+1))]
  colfunc_vec_2 = colfunc(100)[-c((round(100/2)-5):(100/2)-1, ((100/2)+1):((100/2)+1+5))]
  
  svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/DMRs_all_cntx_TEs/legends.svg"), width = 6, height = 5.75, family = "serif")
  #par(mfrow=c(3,1))
  layout(matrix(1:5,ncol=1,nrow=5),height = c(1,0.5,0.15,0.15,0.15))
  
  par(mar = c(0.1, 0.1, 0.1, 0.1) )
  circos.par("track.height" = 0.7, "canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1), "gap.degree" = 0, "clock.wise" = FALSE)
  circos.initialize(factors = as.character(1:tracks_total_size), xlim = c(0, 1)) 
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, area = T, border = T, col = "#FF000080")
  circos.lines(seq(0,1,1/(length(rndm_norm_hypo)-1)), rndm_norm_hypo, area = T, border = T, col = "#304ed180")
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  circos.text(text_tracks_pos, x = 0.1, y = 0.5, labels = "CG context", facing = "downward", ad = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, area = T, border = T, col = "#FF000080")
  circos.lines(seq(0,1,1/(length(rndm_norm_hypo)-1)), rndm_norm_hypo, area = T, border = T, col = "#304ed180")
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  circos.text(text_tracks_pos, x = 0.1, y = 0.5, labels = "CHG context", facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, area = T, border = T, col = "#FF000080")
  circos.lines(seq(0,1,1/(length(rndm_norm_hypo)-1)), rndm_norm_hypo, area = T, border = T, col = "#304ed180")
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  circos.text(text_tracks_pos, x = 0.1, y = 0.5, labels = "CHH context", facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.15) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_TE)-1)), rndm_norm_TE, area = T, border = T, col = "#fcba0330")
  #circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  circos.text(text_tracks_pos, x = 0.1, y = 0.5, labels = "TEs density", facing = "downward", adj = c(0, 0.5))
  
  # DMRs
  legend(-0.09,0.125,
         legend=c(substitute(paste(bold("Hyper"),"-DMRs")),
                  substitute(paste(bold("Hypo"),"-DMRs")),
                  "Shared"),
         fill = c("#FF000095","#304ed195","#5e0d3d99"),
         bty="n")
  
  circos.clear()
  dev.off()
}
