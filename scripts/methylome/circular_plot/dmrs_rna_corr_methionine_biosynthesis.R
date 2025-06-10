library(dplyr)
library(writexl)
library(stringr)
library(RColorBrewer)
library(rtracklayer)
library(circlize)
library(org.At.tair.db)
library(KEGGREST)

#####
# in the plot, each colored line (group line color) is 111.111Kbp (inserts 27 times in 3Mbp)
#####


### GFF and TE files
{
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
  transcript2merge = as.data.frame(transcripts)[,c("seqnames","start","end","ID","gene")]
  names(transcript2merge)[4] = "transcript_id"
  transcript2merge$transcript_id = gsub("transcript:","",transcript2merge$transcript_id)
  
  
  ############# TE file #############
  TE = read.csv("P:/TAIR10.1/TAIR10_Transposable_Elements.txt", sep = "\t")
  TE_4_dens = TE[,c(1,3,4)]
  TE_4_dens$Transposon_Name = gsub(paste0("AT1TE.*"), "NC_003070.9",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT2TE.*"), "NC_003071.7",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT3TE.*"), "NC_003074.8",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT4TE.*"), "NC_003075.7",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT5TE.*"), "NC_003076.8",TE_4_dens$Transposon_Name)
} 

for (treatment in c("mto1","mto3")) {

  # RNAseq res file
  RNA_file_0 = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/met23/",treatment,"_vs_wt/all.transcripts.",treatment,"_vs_wt.DE.csv"))
  RNA_file = RNA_file_0[RNA_file_0$padj < 0.05, c("transcript_id","log2FoldChange")] %>% na.omit()
  #RNA_file = RNA_file_0[, c("transcript_id","log2FoldChange")]
  names(RNA_file)[2] = "RNA_log2FC"
  
  RNA_transcriptsLoc = merge.data.frame(transcript2merge, RNA_file_0, by = "transcript_id")
  
  
  ######### all groups list
  {

    map_ids = "00270"
    group_names = c("Methionine biosynthesis", "EC 2.1.1.-", "Nucleosome (CC)", "Chromatin remodeling")
    #group_color = rainbow(length(group_names))
    group_color = c("#20873c","#FF7F00", "#999999","#984EA3") # group_color = brewer.pal(n=length(group_names)+3, "Set1")[-c(1,2,6)]
    
    ### add from KEGG map ID   
    tairs_groups_list = list()
    for (t_l in 1:length(map_ids)) {
    pathway_info = keggGet(paste0("map",map_ids[t_l]))
    # all tairs for ath map
    tairs = data.frame(locus_tag = as.list(org.At.tairPATH2TAIR)[[map_ids[t_l]]])
    tairs_groups_list[[t_l]] = merge.data.frame(tairs, RNA_file_0, by = "locus_tag")
    }
     
    
    ### add from RNAseq file
    tairs_groups_list[[2]] = RNA_file_0[grep("2\\.1\\.1\\.",RNA_file_0$EC.number),]
    tairs_groups_list[[3]] = RNA_file_0[grep("nucleosome",RNA_file_0$Gene.Ontology..cellular.component.),]
    tairs_groups_list[[4]] = RNA_file_0[c(grep("chromatin remodeling",RNA_file_0$Gene.Ontology..biological.process.),
                                          grep("chromatin remodeling|chromatin-remodeling",tolower(RNA_file_0$Function..CC.))),]
    
    
    names(tairs_groups_list) = group_names
    
    all_interesting_tairs = data.frame(transcript_id=NULL,padj=NULL,col=NULL,group=NULL)
    for (group_l in 1:length(group_names)) {
      n.group = tairs_groups_list[[group_names[group_l]]][, c("transcript_id","padj")] %>% 
        filter(padj < 0.05) %>% mutate(col = group_color[group_l]) %>% mutate(group = group_names[group_l])
      all_interesting_tairs = rbind(all_interesting_tairs,n.group)
    }
    
    all_interesting_tairs = all_interesting_tairs[,-2]
    
    ##### groups data frame for plot
    Chr_group = merge.data.frame(transcript2merge, all_interesting_tairs, by = "transcript_id")[,c("seqnames","start","end","gene","col")]
    Chr_group$tmp = paste(Chr_group$seqnames,Chr_group$start,Chr_group$end, sep = "XX")
    Chr_group = Chr_group[!duplicated(Chr_group$tmp), -grep("tmp", names(Chr_group))]
    
    ##### find overlapping genes to color in black (if not, will overlay in the plot)
    expanded_gr <- GRanges(seqnames = Rle(Chr_group$seqnames),
                           ranges = IRanges(start = Chr_group$start - 35000,
                                            end = Chr_group$end + 35000),
                           strand = "*")
    for (ol.l in 1:nrow(Chr_group)) {
      overlaps <- findOverlaps(expanded_gr, expanded_gr[ol.l])
      overlaps_df = Chr_group[queryHits(overlaps),]
      if (length(unique(overlaps_df$col)) > 1) {
        Chr_group[queryHits(overlaps),]$col = "black"
      }
    }
  }
  
  for (context in c("CG","CHG","CHH")) {
    for (ann in c("genes","promoters")) {
      
      ann.m = ifelse(ann == "genes", "Genes","Promoters")
      ############# DMRs file
      dmrs_file = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/DMRs_",context,"_",treatment,"_vs_wt.csv"))
      dmrs_file$log2FoldChange = log2(dmrs_file$proportion2 / dmrs_file$proportion1)
      dmrs_file = dmrs_file[,c("seqnames","start","end","log2FoldChange")]
      names(dmrs_file)[4] = "meth_log2FC"
      
      ############# methylome and RNAseq file filtered by corr files #############
      corr_file = read.csv(paste0("P:/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/",treatment,"/",context,"/",ann,".corr.",context,".",treatment,".csv"))
      #corr_file = read.csv(paste0("P:/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_TPM/",treatment,"/",context,"/",ann,".corr.",context,".",treatment,".csv"))
      corr_file = corr_file[corr_file$pval < 0.05, c("transcript_id","locus_tag")]
      #corr_file = corr_file[(corr_file$cor>0.8 | corr_file$cor<(-0.8)), c("transcript_id","locus_tag")]
      
      
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
      
      
      color_palette <- colorRampPalette(c("blue","white", "red"))
      length4norm = nrow(Chr_df)+50
      color_palette_vec = color_palette(length4norm)
      color_palette_vec = color_palette_vec[-c((round(length4norm/2)-24):(length4norm/2), ((length4norm/2)+1):((length4norm/2)+1+24))]
      
      scaled_values_meth <- round(((Chr_meth$val - min(Chr_meth$val)) / (max(Chr_meth$val) - min(Chr_meth$val))) * (length(Chr_meth$val) - 1)) + 1
      Chr_meth$col <- color_palette_vec[scaled_values_meth]
      
      scaled_values_RNA <- round(((Chr_RNA$val - min(Chr_RNA$val)) / (max(Chr_RNA$val) - min(Chr_RNA$val))) * (length(Chr_RNA$val) - 1)) + 1
      Chr_RNA$col <- color_palette_vec[scaled_values_RNA]
      
      
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
      
      # TE size color
      #if (length(overlap_ranges_TE_gr) > 0) {
      #  mcols(overlap_ranges_TE_gr)$col <- "black"
      #  mcols(overlap_ranges_TE_gr)$col[width(overlap_ranges_TE_gr) < 500] <- "gray80"
      #  mcols(overlap_ranges_TE_gr)$col[width(overlap_ranges_TE_gr) > 4000] <- "gray20" 
      #}
      overlap_ranges_TE_df = as.data.frame(overlap_ranges_TE_gr)
      
    
      dir.create(paste0("P:/yonatan/methionine/NGS_merged_results/circular_plot_res/DMRs_RNA_corr_biosynthesis/",treatment))
      dir.create(paste0("P:/yonatan/methionine/NGS_merged_results/circular_plot_res/DMRs_RNA_corr_biosynthesis/",treatment,"/",ann.m))
      #####################################
      ############# the plot #############
      svg(paste0("P:/yonatan/methionine/NGS_merged_results/circular_plot_res/DMRs_RNA_corr_biosynthesis/",treatment,"/",ann.m,"/",treatment,"_DMRs_meth_biosynthesis_",context,"_corr_",ann,"_circular_plots.svg"), width = 7.07, height = 7.40, family = "serif")
      
      circos.genomicInitialize(as.data.frame(gff3.trimmed)[,1:3], sector.names = paste0("Chr ",1:5), axis.labels.cex = 0.75, labels.cex = 1.5)
      
      circos.genomicDensity(list(dmrs_file[dmrs_file$meth_log2FC > 0, 1:3],
                                 dmrs_file[dmrs_file$meth_log2FC < 0, 1:3]),
                            bg.col = "#fafcff", bg.border = NA, count_by = "number",
                            col = c("#FF000080","#304ed180"), border = T)
      
      circos.genomicTrackPlotRegion(Chr_group[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, panel.fun = function(region, value, ...) {
        region_colors <- Chr_group$col[which(Chr_group[,1] == get.current.chromosome())]
        circos.genomicLines(region, 1, type = "h", col = region_colors, lwd = 1.5)
      }, track.height = 0.06)
      
      circos.genomicTrackPlotRegion(Chr_meth[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, panel.fun = function(region, value, ...) {
        region_colors <- Chr_meth$col[which(Chr_meth[,1] == get.current.chromosome())]
        circos.genomicLines(region, 1, type = "h", col = region_colors, lwd = 1.5)
      }, track.height = 0.06)
      
      circos.genomicTrackPlotRegion(Chr_RNA[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, panel.fun = function(region, value, ...) {
        region_colors <- Chr_RNA$col[which(Chr_RNA[,1] == get.current.chromosome())]
        circos.genomicLines(region, 1, type = "h", col = region_colors, lwd = 1.5)
      }, track.height = 0.06)
      
      circos.genomicTrackPlotRegion(overlap_ranges_TE_df[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, panel.fun = function(region, value, ...) {
        #region_colors <- overlap_ranges_TE_df$col[which(overlap_ranges_TE_df[,1] == get.current.chromosome())]
        circos.genomicLines(region, 1, type = "h", col = "black", lwd = 1.5)
      }, track.height = 0.05)
      
      
      #circos.genomicDensity(list(TE_4_dens[1:3], RNA_transcriptsLoc[2:4]),
      #                      bg.col = "#fafcff", bg.border = NA, count_by = "number",
      #                      col = c("#fcba0340","#b2b2b290"), border = T)
      circos.genomicDensity(list(RNA_transcriptsLoc[2:4], TE_4_dens[1:3]),
                            bg.col = "#fafcff", bg.border = NA, count_by = "number",
                            col = c("gray","#fcba0320"), border = T)
      
      circos.clear()
      
      dev.off()
      
      cat(">")
    }
  }
}
############################################################################################################
############################################################################################################

# legends
{
  
  counts_norm_hyper = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  counts_norm_hypo = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  rndm_norm_hyper = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0,10))
  rndm_norm_hypo = c(rep(0,10), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
  #rndm_norm_TE = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
  rndm_norm_TE = c(0, rndm_norm_hypo[-c(1:10)], 0)
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
  
  svg(paste0("P:/yonatan/methionine/NGS_merged_results/circular_plot_res/DMRs_RNA_corr_biosynthesis/legends.svg"), width = 6, height = 5.75, family = "serif")
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
  #circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[1], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.06) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(rndm[1:(length(unique(all_interesting_tairs$col))+1)], rep(1,(length(unique(all_interesting_tairs$col))+1)), type = "h", col = c(unique(all_interesting_tairs$col), "black"), lwd = 1.75)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[2], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = as.character(1:tracks_total_size), ylim = c(0, 1), bg.border = NA, track.height = 0.06) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(rndm, rep(1,12), type = "h", col = colfunc_vec_1, lwd = 1.5)
  circos.text(text_tracks_pos, x = 0.35, y = 0.75, labels = legend_titles_circular[3], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.06) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(rndm, rep(1,12), type = "h", col = colfunc_vec_1[12:1], lwd = 1.5)
  circos.text(text_tracks_pos, x = 0.4, y = 0.75, labels = legend_titles_circular[4], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.05) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(rndm[1:4], rep(1,4), type = "h", col = "black", lwd = 1.5)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[5], facing = "downward", adj = c(0, 0.5))
  
  #circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  #circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  #circos.lines(seq(0,1,1/14), rnorm(15,0.5, 0.125), area = T, border = F)
  #circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[6], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_genes)-1)), rndm_norm_genes, area = T, border = T, col = "#b2b2b2")
  circos.lines(seq(0,1,1/(length(rndm_norm_TE)-1)), rndm_norm_TE, area = T, border = T, col = "#fcba0330")
  #circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  #circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[6], facing = "downward", adj = c(0, 0.5))
  
  # DMRs
  legend(0,0.99,
         legend=c(substitute(paste(bold("Hyper"),"-methylated DMRs")),
                  substitute(paste(bold("Hypo"),"-methylated DMRs")),
                  "Shared"),
         fill = c("#FF000095","#304ed195","#5e0d3d99"),
         bty="n")
  #title=expression(bold("DMRs direction")))
  #title="DMRs direction")
  
  # values
  legend_image <- as.raster(matrix(colfunc_vec_2, ncol=1))
  #rasterImage(legend_image, 0.035, 0.605, 0.055,0.685)
  #rect(0.035, 0.605, 0.055,0.685, border="black")
  #text(x=0.04, y = c(0.585,0.705), labels = c(expression(bold(" -1")),expression(bold(" 1"))), cex=0.8)
  rasterImage(legend_image, 0.035, 0.59, 0.055, 0.705)
  rect(0.035, 0.59, 0.055,0.705, border="black")
  text(x=0.0175, y = c(0.605,0.69), labels = c(expression(bold("-1")),expression(bold(" 1"))), cex=0.8)
  
  # genes/TEs
  legend(0,0.45,
         legend=c("sig. Genes","TEs"),
         fill = c("#b2b2b2","#fcba0360"),
         bty="n")
  
  circos.clear()
  
  # group color
  par(mar = c(0, 0.1, 0, 0.1) )
  plot.new()
  legend(0.1,1,#0.2,0.6,
         legend = c(unique(all_interesting_tairs$group),"Overlap"),
         fill = c(unique(all_interesting_tairs$col), "black"),
         bty="n",
         title=expression(bold("Genes Group")),
         title.adj=0, title.cex = 1.1)
  
  # TE size
#  par(mar = c(0, 0.1, 0, 0.1) )
#  plot.new()
#  legend(0.1,1,#0.2,0.6,
#         legend = c("TE","Short TE (<500bp)", "Long TE (>4kbp)"),
#         fill = c("black","gray80","gray20"),
#         bty="n",
#         title=expression(bold("TE size")),
#         title.adj=0, title.cex = 1.1)
  
  dev.off()
  
}
