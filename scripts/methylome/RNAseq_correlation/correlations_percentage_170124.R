for (treatment in c("mto1","mto3")) {
  for (ann in c("genes","promoters")) {
    
    ann.2 = ifelse(ann == "genes", "Genes", "Promoters")
    #tairs = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/description_file_161123.csv")$locus_tag
    tairs_CG = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/",treatment,"_vs_wt/genome_annotation/CG/",ann.2,"_CG_genom_annotations.csv"))$locus_tag
    tairs_CHG = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/",treatment,"_vs_wt/genome_annotation/CHG/",ann.2,"_CHG_genom_annotations.csv"))$locus_tag
    tairs_CHH = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/",treatment,"_vs_wt/genome_annotation/CHH/",ann.2,"_CHH_genom_annotations.csv"))$locus_tag
    
    tairs_CG = length(unique(tairs_CG[grep("^AT[0-5]G", tairs_CG)]))
    tairs_CHG = length(unique(tairs_CHG[grep("^AT[0-5]G", tairs_CHG)]))
    tairs_CHH = length(unique(tairs_CHH[grep("^AT[0-5]G", tairs_CHH)]))
    
    ############# methylome and RNAseq file filtered by corr files #############
    corr_CG = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/",treatment,"/CG/",ann,".corr.CG.",treatment,".csv"))
    corr_CHG = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/",treatment,"/CHG/",ann,".corr.CHG.",treatment,".csv"))
    corr_CHH = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/",treatment,"/CHH/",ann,".corr.CHH.",treatment,".csv"))
    
    corr_CG = corr_CG[corr_CG$pval < 0.05,]
    corr_CHG = corr_CHG[corr_CHG$pval < 0.05,]
    corr_CHH = corr_CHG[corr_CHH$pval < 0.05,]
    
    
    CG_length = c(length(unique(corr_CG[corr_CG$cor < 0,"locus_tag"])),
                  length(unique(corr_CG[corr_CG$cor > 0,"locus_tag"])))
    CHG_length = c(length(unique(corr_CHG[corr_CHG$cor < 0,"locus_tag"])),
                   length(unique(corr_CHG[corr_CHG$cor > 0,"locus_tag"])))
    CHH_length = c(length(unique(corr_CHH[corr_CHH$cor < 0,"locus_tag"])),
                   length(unique(corr_CHH[corr_CHH$cor > 0,"locus_tag"])))
    
    
    CG_pres = (CG_length / tairs_CG) * 100
    CHG_pres = (CHG_length / tairs_CHG) * 100
    CHH_pres = (CHH_length / tairs_CHH) * 100
    
    #### Percentage by DMRs file or sig. genes. not all!!!!
    
    plot_df <- data.frame(CG = CG_pres, CHG = CHG_pres, CHH = CHH_pres)
    rownames(plot_df) <- c("Negative","Positive")
    
    y_lim_perc = max(plot_df) / 9
    y_lim_axis_0 = max(c(sum(plot_df[,1]),sum(plot_df[,2]),sum(plot_df[,3])))
    y_lim_axis = seq(0,95,5)[which(seq(0,95,5) >= y_lim_axis_0)][1]
    
    dir.create("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/corr_percentage_barplots/", showWarnings = F)
    svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/corr_percentage_barplots/",treatment,"_",ann,"_corr_percentage.svg"), width = 2.8, height = 3.85, family = "serif")
    barplot(as.matrix(plot_df), 
            col=c("#6c96d9","#d96c6c"), 
            border="black", 
            space=0.4, 
            font.axis=2, cex.axis = 1,
            ylim = c(0,y_lim_axis),
            xlab="", ylab = "Correlations (%)")
    
    text(0.9,c(plot_df$CG[1]-y_lim_perc, sum(plot_df$CG)-y_lim_perc), CG_length, cex=0.8) 
    text(2.3,c(plot_df$CHG[1]-y_lim_perc, sum(plot_df$CHG)-y_lim_perc), CHG_length, cex=0.8) 
    text(3.7,c(plot_df$CHH[1]-y_lim_perc, sum(plot_df$CHH)-y_lim_perc), CHH_length, cex=0.8) 
    
    
    #text(0.9,c(plot_df$CG[1]/2, plot_df$CG[1]+(plot_df$CG[2]/2)), CG_length, cex=0.8) 
    #text(2.3,c(plot_df$CHG[1]/2, plot_df$CHG[1]+(plot_df$CHG[2]/2)), CHG_length, cex=0.8) 
    #text(3.7,c(plot_df$CHH[1]/2, plot_df$CHH[1]+(plot_df$CHH[2]/2)), CHH_length, cex=0.8) 
    dev.off()
  }
}

### legend ###
svg("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/corr_percentage_barplots/corr_percentage_legend.svg", width = 2.38, height = 3.71, family = "serif")
plot.new()
legend("center", legend = c("Positive","Negative"),
       fill = c("#d96c6c","#6c96d9"),
       bty="n",
       title=expression(bold("Correlation")))
dev.off()
