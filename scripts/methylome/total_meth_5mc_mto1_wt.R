library(ggplot2)

var1 = "wt"
var2 = "mto1"

var1_df = as.data.frame(readxl::read_xlsx("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/old/BSseq_results_191223/total_meth_df.xlsx", sheet = var1))
var2_df = as.data.frame(readxl::read_xlsx("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/old/BSseq_results_191223/total_meth_df.xlsx", sheet = var2))


meth_plot_df = data.frame(type = rep(c("CG", "CHG", "CHH"),2),
                          treatment = c(rep(var1,3), rep(var2,3)),
                          levels = c(var1_df[1:3,2], var2_df[1:3,2]),
                          SD = c(var1_df[1:3,3], var2_df[1:3,3]))
    
    # plot
    level_order = c(var1,var2)
    y_max_plot = max(meth_plot_df$levels)*1.05
    if (nchar(var1) > 6 | nchar(var2) > 6) {leg_horiz=0.75} else {leg_horiz=0.9} # legend position
    g1 <- ggplot(data = meth_plot_df, aes(x = type, y = levels, fill = factor(treatment,level=level_order))) +
      geom_bar(stat = "identity", position = position_dodge(), colour="black") +
      geom_errorbar(aes(ymax = levels + SD, ymin = levels - SD), position = position_dodge(width = 0.9), width = 0.2) +
      scale_fill_manual(values=alpha(c("gray50","#bf6828"),0.8)) +
      theme_classic() +
      theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
            #title = element_text(size = 9, face="bold"),
            axis.title.x = element_blank(),
            #axis.text.x = element_blank(),
            axis.title.y = element_text(size = 12, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12, face="bold"),
            axis.text.y = element_text(size = 10, face="bold"),
            axis.line.x = element_line(linewidth = 1.1),
            axis.line.y = element_line(linewidth = 1.1),
            axis.ticks = element_line(linewidth = 1.1) , 
            axis.ticks.length = unit(0.01, "cm"),
            legend.position = c(leg_horiz, 0.85),
            legend.title = element_blank(),
            legend.text = element_text(size = 9, face = "bold")) + 
      labs(y = "5-mC%") +
      #geom_jitter(color="steelblue4", size=1, alpha=0.9, width = 0.18) +
      scale_y_continuous(expand = c(0,0), limits = c(0,y_max_plot)) 
    
    svg(file = paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/5mC_",var2,"_vs_",var1,".svg"), width = 2.5, height = 4, family = "serif")
    plot(g1)
    dev.off()
