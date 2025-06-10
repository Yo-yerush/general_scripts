library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(ggbreak)

TE = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt", sep = "\t")


TE_DMRs_overlap <- function(TE_family, TE_super_family, TE_df=TE, is.tairs=F) {
  
  TE_df$Transposon_ID = TE_df$Transposon_Name
  TE_df$Transposon_Name = gsub(paste0("AT1TE.*"), "Chr1",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT2TE.*"), "Chr2",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT3TE.*"), "Chr3",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT4TE.*"), "Chr4",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT5TE.*"), "Chr5",TE_df$Transposon_Name)
  
  names(TE_df)[1:4] = c("seqnames","strand","start","end")
  TE_df$strand = ifelse(TE_df$strand == "true","+","-")
  TE_df$strand = "*"
  TE_gr = makeGRangesFromDataFrame(TE_df, keep.extra.columns = T)
  
  if (is.tairs) {
    TE_gr_family = TE_gr[grep(TE_family, TE_gr$Transposon_ID)]
    
  } else if (TE_super_family == "all") { # if its test for super family group compare to all groups
    TE_gr_family = TE_gr[grep(TE_family, TE_gr$Transposon_Super_Family)]
    
  } else {
    TE_gr_family = TE_gr[grep(TE_family, TE_gr$Transposon_Family)]
  }
  
  if (TE_super_family == "all") { # if its test for super family group compare to all groups
    TE_gr_super_family = TE_gr
  } else {
    TE_gr_super_family = TE_gr[TE_gr$Transposon_Super_Family == TE_super_family]
  }
  
  CG = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_191223/mto1_vs_wt/DMRs_CG_mto1_vs_wt.csv")[,c("seqnames","start","end","width","strand","direction")]
  CHG = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_191223/mto1_vs_wt/DMRs_CHG_mto1_vs_wt.csv")[,c("seqnames","start","end","width","strand","direction")]
  CHH = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_191223/mto1_vs_wt/DMRs_CHH_mto1_vs_wt.csv")[,c("seqnames","start","end","width","strand","direction")]
  
  DMRs = rbind(CG,CHG,CHH)
  DMRs$tmp_col = paste(DMRs$seqnames,DMRs$start,DMRs$end, sep = "XXX")
  DMRs = DMRs[!duplicated(DMRs$tmp_col),]
  DMRs = DMRs[,-grep("tmp_col", names(DMRs))]
  DMRs_gr = makeGRangesFromDataFrame(DMRs, keep.extra.columns = T)
  
  # DMRs overlap with family
  family_TE_overlaps <- findOverlaps(TE_gr_family, DMRs_gr)
  family_TE_overlap_ranges <- DMRs_gr[subjectHits(family_TE_overlaps)]
  mcols(family_TE_overlap_ranges) <- cbind(mcols(family_TE_overlap_ranges),mcols(TE_gr_family)[queryHits(family_TE_overlaps),])
  family_TE_overlap_ranges = family_TE_overlap_ranges[!duplicated(family_TE_overlap_ranges$Transposon_ID)]
  
  # DMRs overlap with super family
  super_family_TE_overlaps <- findOverlaps(TE_gr_super_family, DMRs_gr)
  super_family_TE_overlap_ranges <- DMRs_gr[subjectHits(super_family_TE_overlaps)]
  mcols(super_family_TE_overlap_ranges) <- cbind(mcols(super_family_TE_overlap_ranges),mcols(TE_gr_super_family)[queryHits(super_family_TE_overlaps),])
  super_family_TE_overlap_ranges = super_family_TE_overlap_ranges[!duplicated(super_family_TE_overlap_ranges$Transposon_ID)]
  
  
  #### enrichment analysis
  num_family = length(TE_gr_family)
  num_overlap_with_DMRs_in_family = length(family_TE_overlap_ranges)
  
  total_retrotransposons = length(TE_gr_super_family)
  total_overlap_with_DMRs = length(super_family_TE_overlap_ranges)
  
  
  a = num_overlap_with_DMRs_in_family # overlapped TEs in region
  b = num_family - num_overlap_with_DMRs_in_family # *not* overlapped TEs in region
  c = total_overlap_with_DMRs - num_overlap_with_DMRs_in_family # total overlapped TEs (without region)
  d = total_retrotransposons - total_overlap_with_DMRs - b # total *not* overlapped TEs (without region)
  
  contingency_table = matrix(c(a,c,b,d),
                             nrow = 2,
                             dimnames = list(c(TE_family, paste0("Not ",TE_family)),
                                             c("Overlap with DMRs", "Not Overlap with DMRs")))
  
  fisher = fisher.test(contingency_table, alternative = "greater")
  TE_score = as.numeric(fisher$estimate)
  #TE_score = (a/b)/(c/d)
  
  # "gain" or "loss" (direction)  data frame
  direction_family_df = data.frame(hyper = length(family_TE_overlap_ranges[family_TE_overlap_ranges$direction == 1]),
                                   hypo = length(family_TE_overlap_ranges[family_TE_overlap_ranges$direction == -1]))
  direction_super_family_df = data.frame(hyper = length(super_family_TE_overlap_ranges[super_family_TE_overlap_ranges$direction == 1]),
                                         hypo = length(super_family_TE_overlap_ranges[super_family_TE_overlap_ranges$direction == -1]))
  
  # fisher test
  return(list(pValue = fisher$p.value,
              score = TE_score,
              overlapped = num_overlap_with_DMRs_in_family,
              annotated = num_family,
              contingency_table = contingency_table,
              direction_family = direction_family_df,
              direction_super_family = direction_super_family_df))
}


#####################################################


Copia = TE_DMRs_overlap("LTR/Copia", "all")
Gypsy = TE_DMRs_overlap("LTR/Gypsy", "all")
LINE = TE_DMRs_overlap("LINE/L1", "all")
MuDR = TE_DMRs_overlap("DNA/MuDR", "all")
Helitron = TE_DMRs_overlap("RC/Helitron", "all")
SINE = TE_DMRs_overlap("SINE", "all")
Tc1 = TE_DMRs_overlap("DNA/Tc1", "all")

ATCOPIA78 = TE_DMRs_overlap("ATCOPIA78", "LTR/Copia")
ATCOPIA93 = TE_DMRs_overlap("ATCOPIA93", "LTR/Copia")
ATCOPIA51 = TE_DMRs_overlap("ATCOPIA51", "LTR/Copia")
ATCOPIA52 = TE_DMRs_overlap("ATCOPIA52", "LTR/Copia")
atre1 = TE_DMRs_overlap("ATRE1", "LTR/Copia")
romania = TE_DMRs_overlap("ROMANIAT5", "LTR/Copia")
tat1 = TE_DMRs_overlap("TAT1_ATH", "LTR/Gypsy")

ATGP3 = TE_DMRs_overlap("ATGP3", "LTR/Gypsy")

ATHILA = TE_DMRs_overlap("ATHILA", "LTR/Gypsy")
ATHILA2 = TE_DMRs_overlap("ATHILA2", "LTR/Gypsy")
ATHILA4 = TE_DMRs_overlap("ATHILA4", "LTR/Gypsy")
ATHILA5 = TE_DMRs_overlap("ATHILA5", "LTR/Gypsy")
ATHILA6 = TE_DMRs_overlap("ATHILA6", "LTR/Gypsy")
ATHILA6A = TE_DMRs_overlap("ATHILA6A", "LTR/Gypsy")
Athila4_6 = TE_DMRs_overlap("ATHILA4|ATHILA5|ATHILA6", "LTR/Gypsy")

atline = TE_DMRs_overlap("ATLINE", "LINE/L1")
ta1 = TE_DMRs_overlap("TA1", "LINE/L1")

#### ONSEN
ONSEN_df = data.frame(x=paste0("ONSEN",1:8),
                      tair=c("AT1G11265","AT3G61330","AT5G13205","AT1G58140","AT1G48710","AT3G59720","AT1G21945","AT3G32415"),
                      te=c("AT1TE12295","AT3TE92525","AT5TE15240","AT1TE71045","AT1TE59755","AT3TE89830","AT1TE24850","AT3TE54550"))
ONSEN = TE_DMRs_overlap(paste(ONSEN_df$te, collapse = "|"), "LTR/Copia", is.tairs = T)

## EVADE
evd = TE_DMRs_overlap("AT5TE20395", "LTR/Copia", is.tairs = T)

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
