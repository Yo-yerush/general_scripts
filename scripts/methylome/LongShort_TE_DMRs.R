library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(ggbreak)

TE_df = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/TAIR10_Transposable_Elements.txt", sep = "\t")
  
  TE_df$Transposon_ID = TE_df$Transposon_Name
  TE_df$Transposon_Name = gsub(paste0("AT1TE.*"), "NC_003070.9",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT2TE.*"), "NC_003071.7",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT3TE.*"), "NC_003074.8",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT4TE.*"), "NC_003075.7",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT5TE.*"), "NC_003076.8",TE_df$Transposon_Name)
  
  TE_df = TE_df[,-c(5:6)]
  names(TE_df)[1:4] = c("seqnames","strand","start","end")
  TE_df$strand = ifelse(TE_df$strand == "true","+","-")
  TE_df$strand = "*"
  TE_gr = makeGRangesFromDataFrame(TE_df, keep.extra.columns = T) # dont need families characterization hrer
  
  ##########################################

  CG = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/mto1_vs_wt/DMRs_CG_mto1_vs_wt.csv")[,c("seqnames","start","end","width","strand","log2FC")]
  CHG = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/mto1_vs_wt/DMRs_CHG_mto1_vs_wt.csv")[,c("seqnames","start","end","width","strand","log2FC")]
  CHH = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/mto1_vs_wt/DMRs_CHH_mto1_vs_wt.csv")[,c("seqnames","start","end","width","strand","log2FC")]
  
  DMRs_gr = GRangesList(CG = makeGRangesFromDataFrame(CG, keep.extra.columns = T),
                        CHG = makeGRangesFromDataFrame(CHG, keep.extra.columns = T),
                        CHH = makeGRangesFromDataFrame(CHH, keep.extra.columns = T))
  
  ##########################################
  
  # DMRs overlap with TEs
  overlap_fun <- function(te, dmr, is.unique) {
    x <- findOverlaps(te, dmr)
    xx <- te[queryHits(x)]
    mcols(xx) <- cbind(log2FC = mcols(dmr)[subjectHits(x),], size = width(xx), mcols(xx))
    if (is.unique) {
      ## here i only want the average log2FC within TEs. its ok that thee few DMRs on the same TE
      xx = xx[!duplicated(xx$Transposon_ID)]
      ## if its not unique, we can assume that long TEs will overlap more than short (more DMRs on long TEs cous they are long)
    }
    return(xx)
  }

  # short VS long TEs
  range_widths <- width(TE_gr)
  
  long_TEs = TE_gr[range_widths > 4000]
  short_TEs = TE_gr[range_widths < 500]
  
  ##########################################
  
  overlap_TEs = GRangesList(CG = overlap_fun(TE_gr, DMRs_gr[["CG"]], F),
                            CHG = overlap_fun(TE_gr, DMRs_gr[["CHG"]], F),
                            CHH = overlap_fun(TE_gr, DMRs_gr[["CHH"]], F))
  
  overlap_long = GRangesList(CG = overlap_fun(long_TEs, DMRs_gr[["CG"]], F),
                             CHG = overlap_fun(long_TEs, DMRs_gr[["CHG"]], F),
                             CHH = overlap_fun(long_TEs, DMRs_gr[["CHH"]], F))
  
  overlap_short = GRangesList(CG = overlap_fun(short_TEs, DMRs_gr[["CG"]], F),
                              CHG = overlap_fun(short_TEs, DMRs_gr[["CHG"]], F),
                              CHH = overlap_fun(short_TEs, DMRs_gr[["CHH"]], F))
  
  ##########################################
  
  # Calculate the Spearman correlation coefficient
  spearman_correlation <- list(CG = cor(overlap_TEs[["CG"]]$size, overlap_TEs[["CG"]]$log2FC, method = "spearman"),
                               CHG = cor(overlap_TEs[["CHG"]]$size, overlap_TEs[["CHG"]]$log2FC, method = "spearman"),
                               CHH = cor(overlap_TEs[["CHH"]]$size, overlap_TEs[["CHH"]]$log2FC, method = "spearman"))
  
  
  long_spearman_correlation <- list(CG = cor(overlap_long[["CG"]]$size, overlap_long[["CG"]]$log2FC, method = "spearman"),
                                    CHG = cor(overlap_long[["CHG"]]$size, overlap_long[["CHG"]]$log2FC, method = "spearman"),
                                    CHH = cor(overlap_long[["CHH"]]$size, overlap_long[["CHH"]]$log2FC, method = "spearman"))
  
  short_spearman_correlation <- list(CG = cor(overlap_short[["CG"]]$size, overlap_short[["CG"]]$log2FC, method = "spearman"),
                                     CHG = cor(overlap_short[["CHG"]]$size, overlap_short[["CHG"]]$log2FC, method = "spearman"),
                                     CHH = cor(overlap_short[["CHH"]]$size, overlap_short[["CHH"]]$log2FC, method = "spearman"))
  
  # average
  TE_average = list(CG = mean(overlap_TEs[["CG"]]$log2FC),
                    CHG = mean(overlap_TEs[["CHG"]]$log2FC),
                    CHH = mean(overlap_TEs[["CHH"]]$log2FC))
  
  long_average = list(CG = mean(overlap_long[["CG"]]$log2FC),
                      CHG = mean(overlap_long[["CHG"]]$log2FC),
                      CHH = mean(overlap_long[["CHH"]]$log2FC))
  
  short_average = list(CG = mean(overlap_short[["CG"]]$log2FC),
                       CHG = mean(overlap_short[["CHG"]]$log2FC),
                       CHH = mean(overlap_short[["CHH"]]$log2FC))
  
  
  
  ##########################################
  ##########################################
  ##########################################
  ##########################################
  ##########################################
  ##########################################
  ##########################################
  
  
  
  
  
  
  
  
  return(list(pValue = fisher$p.value,
              score = TE_score,
              overlapped = num_overlap_with_DMRs_in_family,
              annotated = num_family,
              contingency_table = contingency_table,
              direction_family = direction_family_df,
              direction_super_family = direction_super_family_df))



  
  
  

df_fun <- function(x) {
  return(data.frame(Score = x$score, pValue = x$pValue, hyper = x$direction_family$hyper, hypo = x$direction_family$hypo))
}

plot_df = rbind(df_fun(ATHILA), df_fun(ATHILA2), df_fun(ATHILA6A), df_fun(ATCOPIA51), df_fun(ATCOPIA78), df_fun(ATGP3), df_fun(ONSEN))
plot_df$family = c("ATHILA","ATHILA2","ATHILA6A","ATCOPIA51","ATCOPIA78","ATGP3","ONSEN")

ggplot(plot_df, aes(x = family, y = Score)) +
  geom_bar(stat = "identity", width = 0.5, fill = "gray60", colour = "black") +
  theme_classic() +
  theme(#panel.spacing = unit(2, "lines"),
    text=element_text(family="serif"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, face="bold"),
    axis.text.y = element_text(face="bold", size=8),
    axis.title.y = element_text(size = 12, face="bold", hjust = 0.65),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 0.75),
    axis.ticks.length.x = unit(-0.2, "cm"),
    axis.ticks.length.y = unit(-0.1, "cm")
    #strip.text = element_blank(),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank()
    #panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  scale_y_break(c(8.25,34), scales = "fixed", ticklabels = c(0,2,4,6,8,35,36)) + 
  labs(y = "Score", x = "") +
  geom_text(aes(label = ifelse(pValue < 0.05, ifelse(pValue < 0.01, ifelse(pValue < 0.001, '***', '**'), '*'), '')),
            vjust = -0.2, size = 4.5) + 
  geom_rect(aes(xmin = 0.5, xmax = nrow(plot_df)+0.5, ymin = 0, ymax = 36),
            fill = "transparent", color = "black", size = 0.75)




#############################################################
############################################################

# to all families from Gypsy and Copia
all_df = data.frame(family=NULL, super_family=NULL, pValue=NULL, Score=NULL, overlapped=NULL, annotated=NULL, direction=NULL)
for (sf in c("LTR/Gypsy","LTR/Copia","LINE/L1")) {
  TE_SF = TE[grep(sf, TE$Transposon_Super_Family),]
  unique_F = unique(TE_SF$Transposon_Family)
  
  for (f in unique_F) {
    loop_run = TE_DMRs_overlap(f, sf)
    all_df = rbind(all_df, data.frame(family=f, super_family=sf,
                                      pValue=loop_run[["pValue"]], Score=loop_run[["score"]],
                                      overlapped=loop_run[["overlapped"]], annotated=loop_run[["annotated"]], 
                                      direction=loop_run[["direction_family"]]))
  }
}


all_df_plot = all_df[order(all_df$Score, decreasing = T),] %>% filter(pValue < 0.05)
row.names(all_df_plot) = all_df_plot$family
all_df_plot$family <- factor(all_df_plot$family, levels = all_df_plot$family)

ggplot(all_df_plot, aes(x = all_df_plot$family, y = Score)) +
  geom_bar(stat = "identity", width = 0.5, fill = "gray60", colour = "black") +
  theme_classic() +
  theme(#panel.spacing = unit(2, "lines"),
    text=element_text(family="serif"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, face="bold"),
    axis.text.y = element_text(face="bold", size=8),
    axis.title.y = element_text(size = 12, face="bold"),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 0.75),
    axis.ticks.length.x = unit(-0.3, "cm"),
    axis.ticks.length.y = unit(-0.05, "cm")
    #strip.text = element_blank(),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  #scale_y_continuous(limits = c(0,7)) +
  labs(y = "Score", x = "") +
  geom_text(aes(label = ifelse(pValue < 0.05, ifelse(pValue < 0.01, ifelse(pValue < 0.001, '***', '**'), '*'), '')),
            vjust = -0.5, size = 4) + 
  geom_rect(aes(xmin = 0.5, xmax = nrow(all_df_plot)+0.5, ymin = 0, ymax = 8),
            fill = "transparent", color = "black", size = 0.75)


#############################################################
############################################################
#############################################################
############################################################
#############################################################
############################################################
#############################################################
############################################################

# percentage of hypermetylated
total_sp = c(Gypsy$annotated,
             Copia$annotated,
             LINE$annotated,
             MuDR$annotated,
             Helitron$annotated)

hyper_sp = c(Gypsy$direction_family$hyper,
             Copia$direction_family$hyper,
             LINE$direction_family$hyper,
             MuDR$direction_family$hyper,
             Helitron$direction_family$hyper)

hypo_sp = c(Gypsy$direction_family$hypo,
            Copia$direction_family$hypo,
            LINE$direction_family$hypo,
            MuDR$direction_family$hypo,
            Helitron$direction_family$hypo)

pValue_sp = c(Gypsy$pValue,
              Copia$pValue,
              LINE$pValue,
              MuDR$pValue,
              Helitron$pValue)

score_sp = c(Gypsy$score,
             Copia$score,
              LINE$score,
              MuDR$score,
              Helitron$score)

level_order =  c("LTR/Gypsy", "LTR/Copia", "LINE/L1", "DNA/MuDR", "RC/Helitron")
plot_df <- data.frame(family = rep(level_order ,2), 
                      direction = c(rep("Hyper",5), rep("Hypo",5)),
                      presentage = c(hyper_sp/total_sp, hypo_sp/total_sp)*100,
                      pValue = c(pValue_sp,rep(1,5)))

#plot_df <- data.frame('LTR/Copia' = c((hyper_copia/total_copia)*100, (hypo_copia/total_copia)*100),
#                      'LTR/Gypsy' = c((hyper_gypsy/total_gypsy)*100, (hypo_gypsy/total_gypsy)*100),
#                      'LINE/L1' = c((hyper_linel1/total_linel1)*100, (hypo_linel1/total_linel1)*100))
#plot_df <- data.frame('opia' = c(hyper_copia, hypo_copia),
#                      'gypsy' = c(hyper_gypsy, hypo_gypsy),
#                      'LINEL1' = c(hyper_linel1, hypo_linel1))
#rownames(plot_df) <- c("Hyper","Hypo")

#y_max_value = max(c(pres_copia, pres_gypsy, pres_linel1))+0.5
y_max_value = (max((hyper_sp+hypo_sp)/total_sp)*100) + 0.5

plot2 = ggplot(plot_df, aes(factor(family, levels = level_order), y = presentage, fill = direction)) +
  geom_bar(stat = "identity", width = 0.5, colour = "black") +
  scale_fill_manual(values = c("Hyper" = "#d96c6c", "Hypo" = "#6c96d9")) +
  theme_classic() +
  theme(#panel.spacing = unit(2, "lines"),
    legend.key.size = unit(0.4, "cm"), 
    text=element_text(family="serif"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, face="bold"),
    axis.text.y = element_text(face="bold", size=8, ),
    axis.title.y = element_text(size = 12, face="bold"),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 0.75),
    axis.ticks.length.x = unit(-0.2, "cm"),
    axis.ticks.length.y = unit(-0.1, "cm")
    
  ) +
  scale_y_continuous(breaks = c(0,5,10,15,20,25)) +
  labs(y = "DMRs Count (%)", x = "", fill = "Direction") +
  geom_text(aes(y = c(((hyper_sp+hypo_sp)/total_sp)*100, rep(0,5)),
                label = ifelse(pValue < 0.05, ifelse(pValue < 0.01, ifelse(pValue < 0.001, '***', '**'), '*'), '')),
            vjust = 0.25, size = 4) + 
  geom_rect(aes(xmin = 0.5, xmax = 5.5, ymin = 0, ymax = 25),#y_max_value),
            fill = "transparent", color = "black", size = 0.75)

svg("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/TEs_SuperFamily_enrichment.svg", width = 2.72, height = 2.91)
plot2
dev.off()
