library(ggplot2)
library(dplyr)

make_group_df <- function(treatment_f, sheet_f, h=3, w) {
  df_x = as.data.frame(readxl::read_xlsx(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - Documents/General/yonatan/methionine/NGS_merged_results/genes_group_results/",treatment_f,"_groups.xlsx"),
                                         sheet = sheet_f)) %>% 
    select(gene_id,Symbol,RNA_log2FC,RNA_pvalue) %>%
    distinct(gene_id, .keep_all = T) %>%
    mutate(Symbol = ifelse(is.na(Symbol), gene_id, Symbol)) %>%
    mutate(Symbol = factor(Symbol, levels = unique(Symbol))) %>%
    mutate(color = ifelse(RNA_log2FC > 0, "red", "blue"),
           significance = case_when(
             RNA_pvalue < 0.001 ~ "***",
             RNA_pvalue < 0.01  ~ "**",
             RNA_pvalue < 0.05  ~ "*",
             TRUE               ~ ""),
           text_pos = ifelse(RNA_log2FC > 0, RNA_log2FC + 0.05, 0.05))
  
  y_max_plot = max(df_x$RNA_log2FC)
  y_min_plot = min(df_x$RNA_log2FC)
  
  adding_max = (y_max_plot + abs(y_min_plot))*0.1
  adding_min = (y_max_plot + abs(y_min_plot))*0.05
  
  g1 <- ggplot(data = df_x, aes(x = reorder(Symbol, -RNA_log2FC), y = RNA_log2FC,  fill = color)) + 
    geom_bar(stat = "identity", position = position_dodge(), colour="black") +
    geom_text(aes(label = significance, y = text_pos), vjust = 0, color = "black") +
    scale_fill_manual(values = c("red" = "#d96c6c", "blue" = "#6c96d9")) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12, face="bold", angle = 90),
      axis.text.y = element_text(size = 10, face="bold"),
      axis.line.x = element_line(linewidth = 1.1),
      axis.line.y = element_line(linewidth = 1.1),
      axis.ticks = element_line(linewidth = 1.1), 
      axis.ticks.length = unit(0.01, "cm"),
      legend.position = 'none'
    ) + 
    labs(y = "Log2(fold change)") +
    scale_y_continuous(expand = c(0, 0), limits = c(y_min_plot - adding_min, y_max_plot + adding_max)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray60")
  
  dir.create(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - Documents/General/yonatan/methionine/rnaseq_23/groups_barPlots/",treatment_f), showWarnings = F)
  svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - Documents/General/yonatan/methionine/rnaseq_23/groups_barPlots/",treatment_f,"/",sheet_f,"_",treatment_f,".svg"), width = w, height = h, family = "serif")
  print(g1)
  dev.off()
}

##############################

for (treatment in c(
  "mto1_vs_wt",
  "mto3_vs_wt",
  "dCGS_vs_EV",
  "SSE_high_vs_EV",
  "SSE_low_vs_EV",
  "SSE_high_vs_SSE_low"
)) {
   make_group_df(treatment, "DNA_methyltransferase", w = 3)
   make_group_df(treatment, "Histone_Lysine_MTs", w = 10)
   make_group_df(treatment, "Royal_Family_Proteins", w = 8)
   make_group_df(treatment, "chromatin_remodeling", w = 17.5)
   make_group_df(treatment, "RdDM_pathway", w = 15)
   make_group_df(treatment, "other_methylation", w = 15)
   make_group_df(treatment, "Cohen_SSE", w = 15)

   cat("\n***", treatment, "***\n")
}

###############################
#
#treatment = "mto1_vs_wt"
#DNA_MT = make_group_df(treatment, "DNA_methyltransferase", h = 2.75, w = 2.65)
#HL_MT = make_group_df(treatment, "Histone_Lysine_MTs", h = 2.75, w = 8.2)
#RF = make_group_df(treatment, "Royal_Family_Proteins", h = 2.75, w = 7.2)
#CHR = make_group_df(treatment, "chromatin_remodeling", h = 2.75, w = 7.5)
#RdDM = make_group_df(treatment, "RdDM_pathway", h = 2.75, w = 2.75)
#other_methlyation = make_group_df(treatment, "other_methylation", h = 2.75, w = 9)
#C_SSE = make_group_df(treatment, "Cohen_SSE", h = 2.75, w = 9)
#
###############################
#
#treatment = "mto3_vs_wt"
#DNA_MT = make_group_df(treatment, "DNA_methyltransferase", h = 2.75, w = 3.5)
#HL_MT = make_group_df(treatment, "Histone_Lysine_MTs", h = 2.75, w = 8.2)
#RF = make_group_df(treatment, "Royal_Family_Proteins", h = 2.75, w = 7.2)
#CHR = make_group_df(treatment, "chromatin_remodeling", h = 2.75, w = 10)
#RdDM = make_group_df(treatment, "RdDM_pathway", h = 2.75, w = 5)
#other_methlyation = make_group_df(treatment, "other_methylation", h = 2.75, w = 9)
#C_SSE = make_group_df(treatment, "Cohen_SSE", h = 2.75, w = 9)
#
###############################
#
#treatment = "dCGS_vs_EV"
#DNA_MT = make_group_df(treatment, "DNA_methyltransferase", h = 2.75, w = 2.65)
#HL_MT = make_group_df(treatment, "Histone_Lysine_MTs", h = 2.75, w = 8.2)
#RF = make_group_df(treatment, "Royal_Family_Proteins", h = 2.75, w = 7.2)
#CHR = make_group_df(treatment, "chromatin_remodeling", h = 2.75, w = 7.5)
#RdDM = make_group_df(treatment, "RdDM_pathway", h = 2.75, w = 2.75)
#other_methlyation = make_group_df(treatment, "other_methylation", h = 2.75, w = 9)
#C_SSE = make_group_df(treatment, "Cohen_SSE", h = 2.75, w = 9)
#
###############################
#
#treatment = "SSE_high_vs_EV"
#DNA_MT = make_group_df(treatment, "DNA_methyltransferase", h = 2.75, w = 2.65)
#HL_MT = make_group_df(treatment, "Histone_Lysine_MTs", h = 2.75, w = 8.2)
#RF = make_group_df(treatment, "Royal_Family_Proteins", h = 2.75, w = 7.2)
#CHR = make_group_df(treatment, "chromatin_remodeling", h = 2.75, w = 7.5)
#RdDM = make_group_df(treatment, "RdDM_pathway", h = 2.75, w = 2.75)
#other_methlyation = make_group_df(treatment, "other_methylation", h = 2.75, w = 9)
#C_SSE = make_group_df(treatment, "Cohen_SSE", h = 2.75, w = 9)
#
###############################
#
#treatment = "SSE_low_vs_EV"
#DNA_MT = make_group_df(treatment, "DNA_methyltransferase", h = 2.75, w = 2.65)
#HL_MT = make_group_df(treatment, "Histone_Lysine_MTs", h = 2.75, w = 8.2)
#RF = make_group_df(treatment, "Royal_Family_Proteins", h = 2.75, w = 7.2)
#CHR = make_group_df(treatment, "chromatin_remodeling", h = 2.75, w = 7.5)
#RdDM = make_group_df(treatment, "RdDM_pathway", h = 2.75, w = 2.75)
#other_methlyation = make_group_df(treatment, "other_methylation", h = 2.75, w = 9)
#C_SSE = make_group_df(treatment, "Cohen_SSE", h = 2.75, w = 9)
#
###############################
#
#treatment = "SSE_high_vs_SSE_low"
#DNA_MT = make_group_df(treatment, "DNA_methyltransferase", h = 2.75, w = 2.65)
#HL_MT = make_group_df(treatment, "Histone_Lysine_MTs", h = 2.75, w = 8.2)
#RF = make_group_df(treatment, "Royal_Family_Proteins", h = 2.75, w = 7.2)
#CHR = make_group_df(treatment, "chromatin_remodeling", h = 2.75, w = 7.5)
#RdDM = make_group_df(treatment, "RdDM_pathway", h = 2.75, w = 2.75)
#other_methlyation = make_group_df(treatment, "other_methylation", h = 2.75, w = 9)
#C_SSE = make_group_df(treatment, "Cohen_SSE", h = 2.75, w = 9)