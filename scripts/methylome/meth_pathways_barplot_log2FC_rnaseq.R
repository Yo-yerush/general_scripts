library(ggplot2)
library(dplyr)

make_group_df <- function(treatment_f, gene_list, gene_list_name, h=3, w) {
  df_x = as.data.frame(readxl::read_xlsx(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_genes.xlsx"),
                                         sheet = treatment_f)) %>% 
    select(gene_id,Symbol,RNA_log2FC,RNA_pvalue) %>%
    distinct(gene_id, .keep_all = T) %>%
    filter(gene_id %in% gene_list) %>%
    mutate(Symbol = ifelse(is.na(Symbol), gene_id, Symbol)) %>%
    #mutate(Symbol = factor(Symbol, levels = unique(Symbol))) %>%
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
  
  g1 <- ggplot(data = df_x, aes(x = reorder(Symbol, Symbol), y = RNA_log2FC, fill = color)) + 
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
    labs(title = gsub("_"," ",treatment_f), y = "Log2(fold change)") +
    scale_y_continuous(expand = c(0, 0), limits = c(y_min_plot - adding_min, y_max_plot + adding_max)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray60")
  
  dir.create("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mutants_selection_220525/barPlots", showWarnings = F)
  svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mutants_selection_220525/barPlots/", gene_list_name, "_",treatment_f,".svg"), width = w, height = h, family = "serif")
  print(g1)
  dev.off()
}

##############################

gene_list_maintain_loop <- c(
"AT1G69770",
"AT5G13960",
"AT3G24360",
"AT4G13460",
"AT5G13960",
"AT5G66750",
"AT3G07610",
"AT3G07660",
"AT3G24360",
"AT4G13460",
"AT3G24340",
"AT4G19020"
)

gene_list_rddm <- c(
"AT2G40030",
"AT4G11130",
"AT2G40040",
"AT5G14620",
"AT5G49160",
"AT3G43920",
"AT2G27040",
"AT2G32940",
"AT5G21150",
"AT1G69770",
"AT2G16390",
"AT3G48670",
"AT2G35160",
"AT5G14920"
)

for (treatment in c(
  "mto1_vs_wt",
  "mto3_vs_wt",
  "dCGS_vs_EV",
  "SSE_high_vs_EV",
  "SSE_low_vs_EV",
  "SSE_high_vs_SSE_low"
)) {
   make_group_df(treatment, gene_list_maintain_loop, "CHG_H3K9me2_loop", w = 2.75)
   make_group_df(treatment, gene_list_rddm, "RdDM_pathway", w = 3.5)

   cat("\n***", treatment, "***\n")
}
