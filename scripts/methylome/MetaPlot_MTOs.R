library(ggplot2)

# function to read and combine files for a given type and group
read_and_combine <- function(group, type, tr_folder, is_TE = is_TE) {
  G_T = ifelse(is_TE, "TEs", "Genes")
  base_path <- paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/old/BSseq_results_040424/results/",tr_folder,"/metaPlots/",G_T,"/metaPlot_tables")
  files <- c("up.stream.csv", "gene.body.csv", "down.stream.csv")
  paths <- paste(base_path, sprintf("%s.%s.%s", group, type, files), sep="/")
  do.call(rbind, lapply(paths, read.csv))
}

#cntx.m = "CHG"
#is_TE = F

var1 = "wt"
var2 = "mto1"
var3 = "mto3"
vars <- c(var1,var2,var3)

v1.col = "gray25"
v2.col = "#bf6828"
v3.col = "#5e5ee6"

for (is_TE in c(T,F)) {
  for (cntx.m in c("CG","CHG","CHH")) {
    
    
    v1.cntx.stream = data.frame(pos = 1:60, Proportion = read_and_combine(var1, cntx.m, "mto1_vs_wt", is_TE)$Proportion)
    v2.cntx.stream = data.frame(pos = 1:60, Proportion = read_and_combine(var2, cntx.m, "mto1_vs_wt", is_TE)$Proportion)
    v3.cntx.stream = data.frame(pos = 1:60, Proportion = read_and_combine(var3, cntx.m, "mto3_vs_wt", is_TE)$Proportion)
    
    
    
    # add column for line color
    v1.cntx.stream$V = "V1"
    v2.cntx.stream$V = "V2"
    v3.cntx.stream$V = "V3"
    
    # add column for line type
    v1.cntx.stream$L = "solid"
    v2.cntx.stream$L = "solid"
    v3.cntx.stream$L = "solid"
    
    v.cntx.stream = rbind(v1.cntx.stream,v2.cntx.stream,v3.cntx.stream)
    
    min_value = min(v.cntx.stream$Proportion)
    max_value = max(v.cntx.stream$Proportion)
    q1_value = min_value+((max_value-min_value)/3)
    q2_value = min_value+((max_value-min_value)/3)+((max_value-min_value)/3)
    #middle_value = round(mean(c(min_value, max_value)), 2)
    q_value = min_value+((max_value-min_value)/4)
    
    #if (is_TE) {
    legend_labels = c(var1, paste0("\n",var2),paste0("\n\n",var3))
    #} else {
    #  legend_labels = c(paste0("\n",var1),paste0("\n\n",var2))
    #}
    
    main_title = ifelse(is_TE, "TEs", "Gene body")
    breaks_and_labels <- list(breaks = c(1, 20, 40, 60), labels = c("    -2kb", "TSS", "TTS", "+2kb    "))
    
    plot_out = ggplot(data = v.cntx.stream, aes(x = pos, y = Proportion, color = V, group = V)) +
      geom_vline(xintercept = c(20, 40), colour = "gray", linetype = "solid", size = 0.5) +
      geom_line(size = 0.5, aes(linetype = L)) +
      
      #scale_color_manual(values = c("V1" = "#e37320", "V2" = "#0072B2")) +
      scale_color_manual(values = c("V1" = v1.col, "V2" = v2.col, "V3" = v3.col)) +
      #scale_linetype_manual(values = c("V1" = "dashed", "V2" = "solid")) +
      theme_classic() +
      labs(title = main_title,
           x = "",
           y = paste0(cntx.m," methylation")) +
      theme(legend.position = "none",
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 9)
      ) +
      scale_x_continuous(breaks = breaks_and_labels$breaks, labels = breaks_and_labels$labels, expand = expand_scale(add = c(0,0))) +
      scale_y_continuous(breaks = c(min_value, q1_value, q2_value, max_value),
                         labels = c(round(min_value, 3), round(q1_value, 3),
                                    round(q2_value, 3), round(max_value, 3))) +
      annotate("text",
               x = 3,
               y = ifelse(is_TE,max_value,q1_value*1.02),
               label = legend_labels,
               hjust = 0, vjust = 0.75, size = 3.25, 
               color = c(v1.col,v2.col,v3.col), fontface = "bold")
    
    
    # if (is_TE) {
    #    svg(paste0("TEs_",cntx.m,"_metaPlot_",var2,"_vs_",var1,".svg"), width = 1.88, height = 1.94, family = "serif")
    #  } else {
    #    svg(paste0("Genes_",cntx.m,"_metaPlot_",var2,"_vs_",var1,".svg"), width = 1.88, height = 1.94, family = "serif")
    #  }
    #par(mar = c(2,2,1,2))
    svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/saved.pics.R/", ifelse(is_TE, "TE","GeneBody"),"_",cntx.m,"_tmp.svg"), width = 2, height = 2.05, family = "serif")
    print(plot_out)
    dev.off()
    
  }
}
