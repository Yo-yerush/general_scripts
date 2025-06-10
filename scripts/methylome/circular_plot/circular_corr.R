library(rtracklayer)
library(dplyr)
library(circlize)

trimm_Chr <- function(gr_obj) {
  remove_seqnames = c("NC_000932.1","NC_037304.1")
  gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
  seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
  return(sort(gr_obj))
}

############# gff3 fNULL############# gff3 file #############
gff3 = import.gff3("P:/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
gff3.trimmed = trimm_Chr(gff3)
#genes = gff3.trimmed[which(gff3.trimmed$type == "gene")]
transcripts = gff3.trimmed[which(gff3.trimmed$type == "transcript")]

############# TE file #############
TE = read.csv("P:/TAIR10.1/TAIR10_Transposable_Elements.txt", sep = "\t")
TE_4_dens = TE[,c(1,3,4)]
TE_4_dens$Transposon_Name = gsub(paste0("AT1TE.*"), "NC_003070.9",TE_4_dens$Transposon_Name)
TE_4_dens$Transposon_Name = gsub(paste0("AT2TE.*"), "NC_003071.7",TE_4_dens$Transposon_Name)
TE_4_dens$Transposon_Name = gsub(paste0("AT3TE.*"), "NC_003074.8",TE_4_dens$Transposon_Name)
TE_4_dens$Transposon_Name = gsub(paste0("AT4TE.*"), "NC_003075.7",TE_4_dens$Transposon_Name)
TE_4_dens$Transposon_Name = gsub(paste0("AT5TE.*"), "NC_003076.8",TE_4_dens$Transposon_Name)


for (treatment in c("mto1","mto3")) {
  for (context in c("CG","CHG","CHH")) {
    for (ann in c("genes","promoters")) {
      
      ann.m = ifelse(ann == "genes", "Genes","Promoters")
      ############# methylome and RNAseq file filtered by corr files #############
      corr_file = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/corr_with_methylations/",treatment,"/",context,"/",ann,".corr.",context,".",treatment,".csv"))
      #corr_file = corr_file[corr_file$pval < 0.05, c("transcript_id","locus_tag")]
      corr_file = corr_file[(corr_file$cor>0.8 | corr_file$cor<(-0.8)), c("transcript_id","locus_tag")]
      
      RNA_file = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/met23/",treatment,"_vs_wt/all.transcripts.",treatment,"_vs_wt.DE.csv"))
      RNA_file = RNA_file[RNA_file$padj < 0.05, c("transcript_id","log2FoldChange")]# %>% na.omit()
      names(RNA_file)[2] = "RNA_log2FC"
      
      transcript2merge = as.data.frame(transcripts)[,c("seqnames","start","end","ID")]
      names(transcript2merge)[4] = "transcript_id"
      transcript2merge$transcript_id = gsub("transcript:","",transcript2merge$transcript_id)
      RNA_transcriptsLoc = merge.data.frame(transcript2merge,RNA_file, by = "transcript_id")
      
      
      RNA_filtered = merge.data.frame(RNA_file, corr_file, by = "transcript_id")
      RNA_filtered = RNA_filtered[,-1]
      
      meth_file = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/genome_annotation/",context,"/",ann.m,"_",context,"_genom_annotations.csv"))
      meth_file$log2FoldChange = log2(meth_file$proportion2 / meth_file$proportion1)
      meth_file = meth_file[,c("locus_tag","seqnames","start","end","log2FoldChange")]
      names(meth_file)[5] = "meth_log2FC"
      
      Chr_df = merge.data.frame(meth_file, RNA_filtered, by = "locus_tag")[,c("seqnames","start","end","meth_log2FC","RNA_log2FC","locus_tag")]
      
      Chr_meth = Chr_df[,1:3]
      Chr_meth$val = -1 + ((Chr_df$meth_log2FC - min(Chr_df$meth_log2FC)) * 2) / (max(Chr_df$meth_log2FC) - min(Chr_df$meth_log2FC))
      
      Chr_RNA = Chr_df[,1:3]
      Chr_RNA$val = -1 + ((Chr_df$RNA_log2FC - min(Chr_df$RNA_log2FC)) * 2) / (max(Chr_df$RNA_log2FC) - min(Chr_df$RNA_log2FC))
      
      # first run 'find_interesting_GOs' script
      
      Chr_name = merge.data.frame(Chr_df, all_interesting_tairs, by = "locus_tag")[,c("seqnames","start","end","col")]
      Chr_name$tmp = paste(Chr_name$seqnames,Chr_name$start,Chr_name$end, sep = "XX")
      Chr_name = Chr_name[!duplicated(Chr_name$tmp), -grep("tmp", names(Chr_name))]
      
      
      #Chr_name = Chr_df[,c(1:3, grep("gene|product", names(Chr_df)))]
      #for (rName in 1:nrow(Chr_name)) {
      #  Chr_name$val[rName] = ifelse(!is.na(Chr_name$gene[rName]), Chr_name$gene[rName], Chr_name$product[rName])
      #}
      #Chr_name = Chr_name[,-grep("gene|product", names(Chr_name))]
      
      
      color_palette <- colorRampPalette(c("blue","gray", "red"))
      scaled_values_meth <- round(((Chr_meth$val - min(Chr_meth$val)) / (max(Chr_meth$val) - min(Chr_meth$val))) * (length(Chr_meth$val) - 1)) + 1
      Chr_meth$col <- color_palette(length(Chr_meth$val))[scaled_values_meth]
      
      scaled_values_RNA <- round(((Chr_RNA$val - min(Chr_RNA$val)) / (max(Chr_RNA$val) - min(Chr_RNA$val))) * (length(Chr_RNA$val) - 1)) + 1
      Chr_RNA$col <- color_palette(length(Chr_RNA$val))[scaled_values_RNA]
      
      
      ####### TE and corr overlay #######
      TE_gr = TE_4_dens
      names(TE_gr) = c("seqnames","start","end")
      TE_gr$strand = "*"
      TE_gr = makeGRangesFromDataFrame(TE_gr)
      
      corr_gr = Chr_df[,1:3]
      corr_gr$strand = "*"
      corr_gr = makeGRangesFromDataFrame(corr_gr)
      
      TE_corr_overlaps <- findOverlaps(corr_gr, TE_gr)
      overlap_ranges_TE_gr <- TE_gr[subjectHits(TE_corr_overlaps)]
      #overlap_ranges_corr_gr <- corr_gr[queryHits(TE_corr_overlaps)]
      
      #####################################
      ############# the plot #############
      svg(paste0("P:/yonatan/methionine/circular_plot_res/corr_TE_exp/",treatment,"_",context,"_corr_",ann,"_circular_plots.svg"), width = 7.07, height = 7.40, family = "serif")
      
      circos.genomicInitialize(as.data.frame(gff3.trimmed)[,1:3], sector.names = paste0("Chr ",1:5), axis.labels.cex = 0.75, labels.cex = 1.5)
      
      circos.genomicTrackPlotRegion(Chr_meth[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, panel.fun = function(region, value, ...) {
        region_colors <- Chr_meth$col[which(Chr_meth[,1] == get.current.chromosome())]
        circos.genomicLines(region, 1, type = "h", col = region_colors, lwd = 1.5)
      }, track.height = 0.08)
      
      circos.genomicTrackPlotRegion(Chr_RNA[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, panel.fun = function(region, value, ...) {
        region_colors <- Chr_RNA$col[which(Chr_RNA[,1] == get.current.chromosome())]
        circos.genomicLines(region, 1, type = "h", col = region_colors, lwd = 1.5)
      }, track.height = 0.08)
      
      circos.genomicTrackPlotRegion(as.data.frame(overlap_ranges_TE_gr)[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, panel.fun = function(region, value, ...) {
        circos.genomicLines(region, 1, type = "h", col = "black", lwd = 1.5)
      }, track.height = 0.06)
      
      circos.genomicTrackPlotRegion(Chr_name[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, panel.fun = function(region, value, ...) {
        region_colors <- Chr_name$col[which(Chr_name[,1] == get.current.chromosome())]
        circos.genomicLines(region, 1, type = "h", col = region_colors, lwd = 1.5)
      }, track.height = 0.08)
      
      circos.genomicDensity(RNA_transcriptsLoc[2:4], bg.col = "#fafcfc", bg.border = NA, count_by = "number", track.height = 0.175)
      circos.genomicDensity(TE_4_dens[1:3], bg.col = "#fafcfc", bg.border = NA, count_by = "number")
      circos.clear()
      
      dev.off()
    }
  }
}


############ legend ############
num.v = c(0.5,0.5)
rndm = runif(12, min = 0, max = 1)
legend_titles_circular = c("Methylation","Expression","Overlapping TEs","Overlapping Terms","sig. Genes","TEs")

svg(paste0("P:/yonatan/methionine/circular_plot_res/corr_TE_exp/circular_legend.svg"), width = 8, height = 5, family = "serif")
par(mar = c(1, 2, 0.1, 0.1) )
circos.par("track.height" = 0.7, "canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1), "gap.degree" = 0, "clock.wise" = FALSE)
circos.initialize(factors = as.character(1:20), xlim = c(0, 1)) 

circos.trackPlotRegion(factors = as.character(1:20), ylim = c(0, 1), bg.border = NA, track.height = 0.08) 
circos.updatePlotRegion(sector.index = "6", bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
circos.lines(rndm, rep(1,12), type = "h", col = rep(c("blue","gray","red"),4), lwd = 1.5)
circos.text(sector.index = "5", x = 0, y = 0.5, labels = legend_titles_circular[1], facing = "inside", adj = c(0, 0.5))

circos.trackPlotRegion(factors = 1:20, ylim = c(0, 1), bg.border = NA, track.height = 0.08) 
circos.updatePlotRegion(sector.index = "6", bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
circos.lines(rndm, rep(1,12), type = "h", col = rep(c("gray","red","blue"),4), lwd = 1.5)
circos.text(sector.index = "5", x = 0, y = 0.5, labels = legend_titles_circular[2], facing = "inside", adj = c(0, 0.5))

circos.trackPlotRegion(factors = 1:20, ylim = c(0, 1), bg.border = NA, track.height = 0.06) 
circos.updatePlotRegion(sector.index = "6", bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
circos.lines(rndm[1:6], rep(1,6), type = "h", col = "black")
circos.text(sector.index = "5", x = 0, y = 0.5, labels = legend_titles_circular[3], facing = "inside", adj = c(0, 0.5))

circos.trackPlotRegion(factors = 1:20, ylim = c(0, 1), bg.border = NA, track.height = 0.08) 
circos.updatePlotRegion(sector.index = "6", bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
circos.lines(rndm[1:8], rep(1,8), type = "h", col = brewer.pal(n=8, "Set1"), lwd = 1.5)
circos.text(sector.index = "5", x = 0, y = 0.5, labels = legend_titles_circular[4], facing = "inside", adj = c(0, 0.5))

circos.trackPlotRegion(factors = 1:20, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
circos.updatePlotRegion(sector.index = "6", bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
circos.lines(seq(0,1,1/11), rep(c(num.v+0.25,num.v),3), area = T, border = F)
circos.text(sector.index = "5", x = 0, y = 0.5, labels = legend_titles_circular[5], facing = "inside", adj = c(0, 0.5))

circos.trackPlotRegion(factors = 1:20, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
circos.updatePlotRegion(sector.index = "6", bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
circos.lines(seq(0,1,1/11), rep(c(num.v,num.v+0.25),3), area = T, border = F)
circos.text(sector.index = "5", x = 0, y = 0.5, labels = legend_titles_circular[6], facing = "inside", adj = c(0, 0.5))

#legend_titles_circular = paste0("yo",1:6)
#for(i in 1:6) {
#  track_height = 0.08 # Adjust based on your track height
#  y_position = 1 - (i-1) * track_height - track_height / 2 # Adjust for the position of the title
#  circos.text(sector.index = "5", x = 0.5, y = y_position, labels = legend_titles_circular[i], facing = "inside", adj = c(0.5, 0.5))
#}

circos.clear()
dev.off()
#############################################################
