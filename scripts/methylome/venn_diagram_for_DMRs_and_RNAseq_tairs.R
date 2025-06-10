venn_di <- function(experiment = "mto.1.3",
                    name1 = "mto1_vs_wt",
                    name2 = "mto3_vs_wt",
                    name3 = "mto1_vs_wt", #rnaseq
                    name4 = "mto3_vs_wt", #rnaseq
                    name5 = NA,
                    DMRcaller.dir = "P:/yonatan/methionine/methylome_23/BSseq_results_161123/",
                    RNAseq.dir = "P:/yonatan/methionine/rnaseq_23/met23/",
                    res.dir = "P:/yonatan/methionine/methylome_23/BSseq_results_161123/",
                    DMR_type,
                    context) {
  
  library(dplyr)
  library(ggplot2)
  library(VennDiagram)
  library("ggVennDiagram")
  library(RColorBrewer)
  
  ###
  # load the data
  
  set1 = read.csv(paste0(DMRcaller.dir,name1,"/genome_annotation/",context,"/",DMR_type,"_",context,"_genom_annotations.csv"))
  set2 = read.csv(paste0(DMRcaller.dir,name2,"/genome_annotation/",context,"/",DMR_type,"_",context,"_genom_annotations.csv"))
  set3 = read.csv(paste0(RNAseq.dir,name3,"/all.transcripts.",name3,".DE.csv"))
  set4 = read.csv(paste0(RNAseq.dir,name4,"/all.transcripts.",name4,".DE.csv"))
  
  set3 = set3[set3$padj < 0.05,]
  set4 = set4[set4$padj < 0.05,]
  
  category.names = c(name1, name2, paste0(name3,"_(RNAseq)"), paste0(name4,"_(RNAseq)"))
  
  tair1 = unique(set1$locus_tag)
  tair2 = unique(set2$locus_tag)
  tair3 = unique(set3$locus_tag)
  tair4 = unique(set4$locus_tag)
  
  x = list(tair1[!is.na(tair1)], tair2[!is.na(tair2)], tair3[!is.na(tair3)], tair4[!is.na(tair4)])
  
  
  ###
  
  # venn diagram position and resolution
  if(is.na(name4) == F) {
    category.position = c(-10, 5, 0, 0)
    resolution = 250
    cex = 0.5
  } else if (is.na(name3) == F) {
    category.position = c(-27, 27, 180)
    resolution = 300
    cex = 0.5
  } else {
    category.position = c(0, 0)
    resolution = 300
    cex = 0.5
  }
  
  #fill = c("#440154ff", '#21908dff', '#fde725ff')
  venn_colors = c("red", "blue","yellow", "green", "purple")
  
  venn.diagram(
    x = x,
    category.names = gsub("_"," ",category.names),
    filename = paste0(res.dir,experiment,"_",context,"_",DMR_type,"_VennDiagram.png"),
    disable.logging = T,
    output = TRUE ,
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = resolution,
    compression = "lzw",
    lwd = 1,
    fill = venn_colors[1:length(x)],
    alpha = rep(0.3, length(x)),
    cex = cex,
    fontfamily = "sans",
    cat.cex = 0.4,
    cat.default.pos = "outer",
    cat.pos = category.position,
    cat.fontface = 2,
    cat.fontfamily = "sans"
    #    cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
    #    col=venn_colors,
    #    rotation = 1
  )
  
  
  ###
  
  sets_merge = merge.data.frame(set1,set2, by = "locus_tag")
  
  
  if(is.na(name4) == F) {
    sets_merge.2 = merge.data.frame(sets_merge, set3, by = "locus_tag")
    merge.file = merge.data.frame(sets_merge.2, set4, by = "locus_tag")
  } else if (is.na(name3) == F) {
    merge.file = merge.data.frame(sets_merge, set3, by = "locus_tag")
  } else {
    merge.file = sets_merge
  }
  
  write.csv(merge.file[!duplicated(merge.file$locus_tag),], paste0(res.dir,experiment,"_",context,"_",DMR_type,"_VennDiagram.csv"), row.names=FALSE)
}

for (ann in c("Genes","Promoters")) {
  for (cntx in c("CG","CHG","CHH")) {
    venn_di(DMR_type = ann, context = cntx)
  }
}

