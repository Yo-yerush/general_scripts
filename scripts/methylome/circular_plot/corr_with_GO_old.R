library(dplyr)
library(writexl)
library(stringr)
library(RColorBrewer)
library(rtracklayer)
library(circlize)

# part1 = find interesting GOs
{
  yo=read.csv("P:/yonatan/methionine/rnaseq_23/description_file_161123.csv")
  yoGO = yo[,grep("locus_tag|Gene.Ontology",names(yo))]
  
  if (F) {
    GO2Tair = rbind(data.frame(locus_tag = yoGO$locus_tag ,term = yoGO$Gene.Ontology..biological.process.),
                    data.frame(locus_tag = yoGO$locus_tag ,term = yoGO$Gene.Ontology..molecular.function.),
                    data.frame(locus_tag = yoGO$locus_tag ,term = yoGO$Gene.Ontology..cellular.component.))
    GO2Tair = GO2Tair[!GO2Tair$term == "",] %>% na.omit()
    GO2Tair = as.data.frame(do.call(rbind, apply(GO2Tair, 1, function(x) {
      do.call(expand.grid, strsplit(x, "; "))
    })))
    GO2Tair$locus_tag = as.character(GO2Tair$locus_tag)
    GO2Tair$term = as.character(GO2Tair$term)
    GO2Tair$tmp = paste(GO2Tair$locus_tag, GO2Tair$term, sep = " XXX ")
    GO2Tair = GO2Tair[!duplicated(GO2Tair$tmp),-3]
    GO2Tair$term = str_extract(GO2Tair$term, "GO:\\d+")
    write.csv(GO2Tair, "P:/TAIR10.1/GO_2_Tair.csv", row.names = F)
  }
  GO2Tair = read.csv("P:/TAIR10.1/GO_2_Tair.csv")
  
  
  bp_df = data.frame(x=yoGO$Gene.Ontology..biological.process.) %>% na.omit()
  bp_df = as.data.frame(do.call(rbind, apply(bp_df, 1, function(x) {
    do.call(expand.grid, strsplit(x, "; "))
  })))
  bp = as.character(unique(bp_df$x))
  
  mf_df = data.frame(x=yoGO$Gene.Ontology..molecular.function.) %>% na.omit()
  mf_df = as.data.frame(do.call(rbind, apply(mf_df, 1, function(x) {
    do.call(expand.grid, strsplit(x, "; "))
  })))
  mf = as.character(unique(mf_df$x))
  
  cc_df = data.frame(x=yoGO$Gene.Ontology..cellular.component.) %>% na.omit()
  cc_df = as.data.frame(do.call(rbind, apply(cc_df, 1, function(x) {
    do.call(expand.grid, strsplit(x, "; "))
  })))
  cc = as.character(unique(cc_df$x))
  
  GO_df = rbind(data.frame(type = rep("biological_process",length(bp)),term = bp),
                data.frame(type = rep("molecular_function",length(mf)),term = mf),
                data.frame(type = rep("cellular_component",length(cc)),term = cc))
  
  all_methylation = GO_df[grep("methylation",GO_df$term, ignore.case = TRUE),]
  all_methyltransferase = GO_df[grep("methyltransferase",GO_df$term, ignore.case = TRUE),]
  all_histone = GO_df[grep("histone",GO_df$term, ignore.case = TRUE),]
  all_DNAc5 = GO_df[grep("DNA \\(cytosine-5", GO_df$term, ignore.case = TRUE),]
  all_chromatin = GO_df[grep("chromatin",GO_df$term, ignore.case = TRUE),]
  all_CellWall = GO_df[grep("cell wall",GO_df$term, ignore.case = TRUE),]
  all_methionine = GO_df[grep("methionine",GO_df$term, ignore.case = TRUE),]
  all_others = GO_df[grep("epigenetic|methylated|methyl-Cp|COMPASS",GO_df$term, ignore.case = TRUE),]
  
  all_DFs = list(methyltransferase = all_methyltransferase,
                 DNA_5mC = all_DNAc5,
                 histone = all_histone,
                 chromatin = all_chromatin,
                 methylation = all_methylation,
                 cell_wall = all_CellWall,
                 methionine = all_methionine,
                 others = all_others)
  #write_xlsx(all_DFs, "P:/yonatan/methionine/interesting_GO_terms/interesting_GO_terms.xlsx")
  
  # extract the go ids for the list
  for (i in 1:length(all_DFs)) {
    all_DFs[[i]]$term = str_extract(all_DFs[[i]]$term, "GO:\\d+")
    for (ii in 1:nrow(all_DFs[[i]])) {
      all_DFs[[i]]$locus_tag[ii] = paste(GO2Tair[grep(all_DFs[[i]]$term[ii], GO2Tair$term), "locus_tag"],
                                         collapse = "; ")
    }
    all_DFs[[i]] = all_DFs[[i]][,-1]
    all_DFs[[i]] = as.data.frame(do.call(rbind, apply(all_DFs[[i]], 1, function(x) {
      do.call(expand.grid, strsplit(x, "; "))
    })))
    all_DFs[[i]]$term = as.character(all_DFs[[i]]$term)
    all_DFs[[i]]$locus_tag = as.character(all_DFs[[i]]$locus_tag)
    all_DFs[[i]]$col = brewer.pal(n=length(all_DFs)+1, "Set1")[-6][i]
  }
  
  all_interesting_tairs = do.call(rbind, all_DFs)[,-1]
  row.names(all_interesting_tairs) = 1:nrow(all_interesting_tairs)
}


# part2 - circular corr plot
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
        
        
        Chr_name = merge.data.frame(Chr_df, all_interesting_tairs, by = "locus_tag")[,c("seqnames","start","end","col")]
        Chr_name$tmp = paste(Chr_name$seqnames,Chr_name$start,Chr_name$end, sep = "XX")
        Chr_name = Chr_name[!duplicated(Chr_name$tmp), -grep("tmp", names(Chr_name))]
        
        
        #Chr_name = Chr_df[,c(1:3, grep("gene|product", names(Chr_df)))]
        #for (rName in 1:nrow(Chr_name)) {
        #  Chr_name$val[rName] = ifelse(!is.na(Chr_name$gene[rName]), Chr_name$gene[rName], Chr_name$product[rName])
        #}
        #Chr_name = Chr_name[,-grep("gene|product", names(Chr_name))]
        
        
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
}


# legends
{
  num.v = c(0.5,0.5)
  rndm = runif(12, min = 0, max = 1)
  legend_titles_circular = c("Methylation","Expression","Overlapping TEs","Overlapping Terms","sig. Genes","TEs")
  tracks_total_size = 12
  legend_tracks_pos = as.character((tracks_total_size/4)+1)
  text_tracks_pos = as.character((tracks_total_size/4))
  
  colfunc <- colorRampPalette(c("red", "white", "blue"))
  colfunc_vec_1 = colfunc(16)[-c((round(16/2)-1):(16/2), ((16/2)+1):((16/2)+1+1))]
  colfunc_vec_2 = colfunc(100)[-c((round(100/2)-5):(100/2)-1, ((100/2)+1):((100/2)+1+5))]
  
  svg(paste0("P:/yonatan/methionine/circular_plot_res/corr_TE_exp/legends.svg"), width = 6, height = 5.75, family = "serif")
  #par(mfrow=c(3,1))
  layout(matrix(1:3,ncol=1,nrow=3),height = c(1,0.25,0.75))
  
  par(mar = c(0.1, 0.1, 0.1, 0.1) )
  circos.par("track.height" = 0.7, "canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1), "gap.degree" = 0, "clock.wise" = FALSE)
  circos.initialize(factors = as.character(1:tracks_total_size), xlim = c(0, 1)) 
  
  circos.trackPlotRegion(factors = as.character(1:tracks_total_size), ylim = c(0, 1), bg.border = NA, track.height = 0.08) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(rndm, rep(1,12), type = "h", col = colfunc_vec_1, lwd = 1.5)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[1], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.08) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(rndm, rep(1,12), type = "h", col = colfunc_vec_1[12:1], lwd = 1.5)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[2], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.06) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(rndm[1:4], rep(1,4), type = "h", col = "black")
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[3], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.08) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(rndm[1:8], rep(1,8), type = "h", col = brewer.pal(n=9, "Set1")[-6], lwd = 1.5)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[4], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/11), rep(c(num.v+0.25,num.v),3), area = T, border = F)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[5], facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/11), rep(c(num.v,num.v+0.25),3), area = T, border = F)
  circos.text(text_tracks_pos, x = 0, y = 0.5, labels = legend_titles_circular[6], facing = "downward", adj = c(0, 0.5))
  
  circos.clear()
  
  #plot.new()
  legend_image <- as.raster(matrix(colfunc_vec_2, ncol=1))
  par(mar = c(0, 10.75, 1, 33) )
  plot(c(0.5,2),c(0,1.1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
  title(main = expression(bold("Values")), line = c(0.5,1), cex=0.9)
  text(x=1.75, y = c(0.05,0.5,0.95), labels = c(-1,0,1), cex=1)
  rasterImage(legend_image, 0, 0, 1,1)
  #box(lty = 1, col = 'black')
  
  par(mar = c(0, 0.1, 0, 0.1) )
  plot.new()
  legend(0.2,0.6,
         legend=names(all_DFs),
         fill = brewer.pal(n=length(all_DFs)+1, "Set1")[-6],
         bty="n",
         title=expression(bold("Terms Group")),
         title.adj=0, title.cex = 1.1)
  #plot.new()
  #  legend(0.2,1.05,
  #         legend=c(1,0,-1),
  #         fill = colorRampPalette(c("red","white", "blue"))(3),
  #        bty="n",
  #         title=expression(bold("Values")),
  #         title.adj=1)
  dev.off()
  
}