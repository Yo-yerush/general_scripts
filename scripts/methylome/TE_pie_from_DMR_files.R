for (context in c("CG","CHG","CHH")) {
  DMRsReplicates_TE_file.0 = paste0("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/mto1_vs_wt/genome_annotation/",context,"/Transposable_Elements_",context,"_genom_annotations.csv")
  DMRsReplicates_TE_file = read.csv(DMRsReplicates_TE_file.0)
  
  # TE pie plot
  TE_Freq = as.data.frame(table(DMRsReplicates_TE_file$Transposon_Super_Family))
  color_vec = c(brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 9, name = "Set1"),"#F7FBFF")
  color_vec = color_vec[1:nrow(TE_Freq)]
  
  svg(file = paste0("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/mto1_vs_wt/genome_annotation/",context,"/Transposon.Super.Family_",context,".svg"), width = 8, height = 5.5, family = "serif")
  par(mar = c(1,4,2,4))
  par(fig=c(0,6,0,10)/10)
  pie(TE_Freq[,2], as.character(TE_Freq[,1]), main = paste0("Transposon Super Family:\nDMRs in ",context," context"),
      col = color_vec)
  par(new=T)
  
  par(mar = c(1,4,2,4))
  par(fig=c(6,10,0,10)/10)
  legend("top", legend = as.character(TE_Freq[,1]), fill = color_vec, bty = "n")
  
  dev.off()
}