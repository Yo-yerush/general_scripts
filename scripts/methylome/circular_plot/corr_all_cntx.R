library(dplyr)
library(RColorBrewer)
library(rtracklayer)
library(circlize)

treatment = "mto1"

### GFF RNAseq, and TEs files
{
  trimm_Chr <- function(gr_obj) {
    remove_seqnames = c("ChrM","ChrC")
    gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
    seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
    return(sort(gr_obj))
  }
  
  ############# gff3 file #############
  gff3 = import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 gff3/TAIR10_GFF3_genes.gff")
  gff3.trimmed = trimm_Chr(gff3)
  max_pos_chr5 = gff3.trimmed[which(gff3.trimmed$type == "chromosome")]@ranges@width[5]
  #genes = gff3.trimmed[which(gff3.trimmed$type == "gene")]
  
  #############  RNAseq res file   #############
  RNA_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/",treatment,"_vs_wt/all.transcripts.",treatment,"_vs_wt.DE.csv")) %>%
    filter(padj < 0.05) %>%
    select(locus_tag, log2FoldChange)
  transcripts = gff3.trimmed[which(gff3.trimmed$type == "mRNA")] %>%
    as.data.frame() %>%
    mutate(locus_tag = as.character(Parent)) %>%
    select(seqnames,start,end,locus_tag) %>%
    distinct(locus_tag, .keep_all = T)
  RNA_transcriptsLoc = merge.data.frame(transcripts, RNA_file, by = "locus_tag")[,c(2:4,1)]
  
  ############# TE file #############
  TE = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt", sep = "\t")
  TE_4_dens = TE[,c("Transposon_Name","Transposon_min_Start","Transposon_max_End","Transposon_Super_Family")]
  TE_4_dens$Transposon_ID = TE_4_dens$Transposon_Name
  for (chr.i in 1:5) {
    TE_4_dens$Transposon_Name = gsub(paste0("AT",chr.i,"TE.*"), paste0("Chr",chr.i),TE_4_dens$Transposon_Name)
  }
  
  ############# TAIR9 gaps #############
  tair9_gaps = rtracklayer::import.gff("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR9/tair9_Assembly_gaps.gff")
  tair9_gaps = resize(tair9_gaps, width = width(tair9_gaps) + 100000, fix = "center") %>%
    reduce() %>%  as.data.frame() %>% .[,1:3]
  
  ############# Heterochromatin positions #############
  heteroChr = data.frame(Chr = paste0("Chr",c(1:5)),
                         start = c(12500000, 1250000, 11000000, 1666667, 9444444),
                         end = c(17500000, 7500000, 16250000, 7000000, 15000000))
  
  ############# Centromere positions #############
  cenChr = data.frame(Chr = paste0("Chr",c(1:5)),
                      start = c(14476796, 3462971, 13780083, 3177188, 11207348),
                      end = c(15081019, 3650512, 14388500, 3248799, 11555278))
  #cenChr = data.frame(Chr = paste0("Chr",c(1:3,3,3:5)),
  #                    start = c(15086046, 3607930, 13799418, 13587787, 14208953, 3956022, 11725025),
  #                    end = c(15087045, 3608929, 13800417, 13588786, 14209952, 3957021, 11726024))
  #cenChr = data.frame(Chr = paste0("Chr",c(1:5)),
  #                    start = c(15086046, 3607930, 13799418, 3956022, 11725025),
  #                    end = c(15087045, 3608929, 13800417, 3957021, 11726024))
  #cenChr$start = cenChr$start - 4.9e4
  #cenChr$end = cenChr$end + 4.9e4
  
  ###############
  ### Calculate densities for RNAseq and TEs
  density_RNA <- genomicDensity(RNA_transcriptsLoc[,1:3], window.size = 1e6, count_by = "number")
  density_TE <- genomicDensity(TE_4_dens[,1:3], window.size = 1e6, count_by = "number")
  # new scaling
  density_RNA$value = density_RNA$value / max(density_RNA$value) * 0.4
  density_TE$value = density_TE$value / max(density_TE$value)
}


for (ann in c("genes","promoters")) {
  
  ann.m = ifelse(ann == "genes", "Genes","Promoters")
  
  
  corr_df_fun <- function(context) {
    
    ############# methylome and RNAseq file filtered by corr files #############
    corr_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/",treatment,"/",context,"/",ann,".corr.",context,".",treatment,".csv"))
    #corr_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_TPM/",treatment,"/",context,"/",ann,".corr.",context,".",treatment,".csv"))
    corr_file = corr_file %>% filter(pval < 0.05) %>% select(-padj)
    #corr_file = data.frame(locus_tag = corr_file[corr_file$pval < 0.05,])
    #corr_file = corr_file[(corr_file$cor>0.8 | corr_file$cor<(-0.8)), c("transcript_id","locus_tag")]
    
    
    RNA_filtered = merge.data.frame(RNA_file, corr_file, by = "locus_tag")
    names(RNA_filtered) = c("gene_id","RNA_log2FC","cor","cor_pval")
    
    meth_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/",treatment,"_vs_wt/genome_annotation/",context,"/",ann.m,"_",context,"_genom_annotations.csv"))
    #meth_file = meth_file[,c("gene_id","seqnames","start","end","log2FC","Symbol")]
    #names(meth_file)[5] = "meth_log2FC"
    meth_file = meth_file %>% dplyr::rename(meth_log2FC = "log2FC")
    
    final_diltered_df_0 = merge.data.frame(meth_file, RNA_filtered, by = "gene_id")
    # final data frame for save
    final_diltered_df = final_diltered_df_0[,-c(2:7,grep("sumReads|proportion|cytosinesCount|pValue|direction|regionType",
                                                         names(final_diltered_df_0)))] %>%
      relocate(RNA_log2FC, cor, cor_pval, .after = meth_log2FC)
    
    # data frames for plots
    Chr_df = final_diltered_df_0[,c("seqnames","start","end","meth_log2FC","RNA_log2FC","gene_id")]
    
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
    
    return(list(meth = Chr_meth,
                RNA = Chr_RNA,
                df = final_diltered_df))
  }
  
  
  dir.create(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/",treatment,"_paper_DMRs_RNA_corr"), showWarnings = F)
  dir.create(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/",treatment,"_paper_DMRs_RNA_corr/",treatment), showWarnings = F)
  dir.create(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/",treatment,"_paper_DMRs_RNA_corr/",treatment,"/",ann.m), showWarnings = F)
  #####################################
  ############# the plot #############
  #svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/mto1_paper_DMRs_RNA_corr/",treatment,"/",ann.m,"/",treatment,"_DMRs_DEGs_corr_",ann,"_circular_plots.svg"),
  #width = 4.25, height = 4.25, family = "serif")
  
  circos.par(gap.degree = c(rep(4,4),35), start.degree = 90)
  circos.genomicInitialize(as.data.frame(gff3.trimmed)[,1:3], sector.names = paste0("Chr ",1:5), axis.labels.cex = 0.4, labels.cex = 1.35)
  
  for (cntx.i in c("CG","CHG","CHH")) {
    
    Chr_meth = corr_df_fun(cntx.i)$meth
    Chr_RNA = corr_df_fun(cntx.i)$RNA
    
    # save the results data frame
    write.csv(corr_df_fun(cntx.i)$df, 
              paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/mto1_paper_DMRs_RNA_corr/",
                     treatment,"/",ann.m,"/",treatment,"_DMRs_DEGs_corr_",ann,"_",cntx.i,".csv"),
              row.names = F)
    
    # DMRs
    circos.genomicTrackPlotRegion(Chr_meth[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, bg.lty = 2, panel.fun = function(region, value, ...) {
      region_colors <- Chr_meth$col[which(Chr_meth[,1] == get.current.chromosome())]
      circos.genomicLines(region, 1, type = "h", col = region_colors, lwd = 1.15)
      ## draw top border
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim") * 1.08 # adjust ylim for top line border
      circos.rect(xlim[1], ylim[2], xlim[2], ylim[2], border = "gray", lwd = 0.85, lty = 1)
      ##
    }, track.height = 0.075, track.margin=c(0, 0))
    
    # DEGs
    circos.genomicTrackPlotRegion(Chr_RNA[,1:3], ylim = c(0,1), bg.col = "#fafcfc", bg.border = NA, bg.lty = 2, panel.fun = function(region, value, ...) {
      region_colors <- Chr_RNA$col[which(Chr_RNA[,1] == get.current.chromosome())]
      circos.genomicLines(region, 1, type = "h", col = region_colors, lwd = 1.15)
      ## draw bottom border
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      ylim = -(ylim[2] * 0.08) # adjust ylim for bottom line border
      circos.rect(xlim[1], ylim[1], xlim[2], ylim[1], border = "gray", lwd = 0.85, lty = 1)
      ##
    }, track.height = 0.075, track.margin=c(0.03, 0))
    
    ### y-axis labels
    circos.text("Chr1", x = 0, y = 0.75, labels = paste0(cntx.i," "), facing = "downward", cex = 0.85, adj = c(1,0))
  }
  
  
  # DEGs
  circos.genomicTrackPlotRegion(density_RNA, ylim = c(0,1), bg.border = NA, track.height = 0.2, track.margin = c(0, 0), panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, col = "gray80", border = TRUE, type = "l", area = T)})
  # TEs
  circos.genomicTrackPlotRegion(density_TE, bg.border = NA, track.height = 0.05, track.margin = c(0, 0), panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, col = "#fcba0320", border = TRUE, type = "l", area = T,
                        track.index = get.cell.meta.data("track.index")-1 )})
  #circos.genomicDensity(list(RNA_transcriptsLoc[,1:3], TE_4_dens[,1:3]), bg.col = "white", # "#fafcff",
  #                      bg.border = NA, count_by = c("percent", "number"), ylim.force = T,
  #                      col =c("gray80","#fcba0320"), border = T, track.height = 0.2, track.margin=c(0, 0))
  circos.clear()
  dev.off()
  
  cat(">")
}



####################################################
# legends
if (F) {
  colfunc <- colorRampPalette(c("red", "white", "blue"))
  colfunc_vec_1 = colfunc(16)[-c((round(16/2)-1):(16/2), ((16/2)+1):((16/2)+1+1))]
  colfunc_vec_2 = colfunc(100)[-c((round(100/2)-5):(100/2)-1, ((100/2)+1):((100/2)+1+5))]
  
  #counts_norm_hyper = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  #counts_norm_hypo = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  #rndm_norm_hyper = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0,10))
  #rndm_norm_hypo = c(rep(0,10), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
  #
  ##rndm_norm_TE = c(rep(0,10),(counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0,10))
  #rndm_norm_TE = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
  ##rndm_norm_TE = c(rep(0,10), rndm_norm_hypo[-c(1:10, (length(rndm_norm_hypo)-10):length(rndm_norm_hypo))], rep(0,10))
  #
  rndm_norm_TE = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  rndm_norm_TE = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)))
  rndm_norm_genes = (1-rndm_norm_TE)/1.4
  
  n.rndm = 8
  rndm = runif(n.rndm, min = 0, max = 1)
  legend_titles_circular = c("DMRs","  Group sig. genes","Methylation values","Expression values","  Overlapping TEs","Genome annotation")
  tracks_total_size = 24
  legend_tracks_pos = as.character((tracks_total_size/4)+1)
  text_tracks_pos = as.character((tracks_total_size/4))
  svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/mto1_paper_DMRs_RNA_corr/legend_corr_circular_plots.svg"),
      width = 6, height = 5.75, family = "serif")
  
  circos.par("track.height" = 0.7, "canvas.xlim" = c(-0.5, 1), "canvas.ylim" = c(0, 1), "gap.degree" = 0, "clock.wise" = FALSE)
  circos.initialize(factors = as.character(1:tracks_total_size), xlim = c(0, 1)) 
  
  circos.trackPlotRegion(factors = as.character(1:tracks_total_size), ylim = c(0, 1), bg.border = NA, track.height = 0.06) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(rndm, rep(1,n.rndm), type = "h", col = colfunc_vec_1, lwd = 1.5)
  circos.text(text_tracks_pos, x = 0.35, y = 0.75, labels = legend_titles_circular[3], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.06) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(rndm, rep(1,n.rndm), type = "h", col = colfunc_vec_1[12:1], lwd = 1.5)
  circos.text(text_tracks_pos, x = 0.4, y = 0.75, labels = legend_titles_circular[4], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_genes)-1)), rndm_norm_genes, area = T, border = T, col = "#b2b2b2")
  circos.lines(seq(0,1,1/(length(rndm_norm_TE)-1)), rndm_norm_TE, area = T, border = T, col = "#fcba0330")
  #circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  #circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[6], facing = "downward", adj = c(0, 0.5))
  
  # values
  #par(mar = rep(0,4))
  #plot(1, type="n", axes=FALSE, xlab="", ylab="", main="\nLog2FC", cex.main=1)
  legend_image <- as.raster(matrix(colfunc_vec_2, ncol=1)) 
  rasterImage(legend_image, 0.04, 0.86, 0.06,0.98) # xleft, ybottom, xright, ytop
  rect(0.04, 0.86, 0.06,0.98, border="black") # box for scale
  rect(-0.05, 0.85, 0.099, 1.025, border="black") # box for the small legend
  text(x=0.035, y = 1, labels = expression(bold("log2FC")), cex=0.75)
  text(x=0.0185, y = c(0.88,0.965), labels = c(expression(bold("-1")),expression(bold(" 1"))), cex=rep(0.75,2))
  
  # genes/TEs
  legend(0,0.825,
         legend=c("sig. Genes","TEs"),
         fill = c("#b2b2b2","#fcba0360"),
         bty="n")
  
  circos.clear()
  dev.off()
}
