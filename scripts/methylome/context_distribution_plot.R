library(ggplot2)
library(RColorBrewer)
library(scales)
library(cowplot)

for (treatment in c("mto1","mto3","dCGS")) {
  
  control = ifelse(treatment == "dCGS", "EV", "wt")
  
  CX_distribution_path = "C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/context_distribution_CX/"
  
  color_vec = c("green4", brewer.pal(n=6,name='Blues')[6:3], brewer.pal(n=9,name='YlOrRd')[9:1],"white")
  ann_list = c("Genes","Promoters","CDS","Introns","fiveUTRs","threeUTRs",
               "Transposable_Elements","TEG","Retrotransposons")
  context_vec = c("CG",
                  "CHG","CAG", "CTG", "CCG",
                  "CHH","CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC")
  
  ##### for plot saving
  pdf(paste0(CX_distribution_path,treatment,"_CX_annotations.pdf"),
      width = 20, height = 2.66, family = "serif")
  #####
  
  for (ann.l in 1:length(ann_list)) {
    
    # Initialize a list to store all plots
    plot_list = list()
    # first plot for 'yaxis' title
    plot_list[[1]] = ggplot(data.frame(x = 1, y = 1), aes(x, y)) + 
      geom_text(aes(label = "5mC Sites"), size = 6, angle = 90, vjust = 0.75) + # Set the text and size
      theme_void() +
      theme(plot.margin = margin(t = 30, r = 0, b = 10, l = 0, unit = "pt"))
    
    for (context_n in 1:length(context_vec)) {
      
      cntx_counts = read.csv(paste0(CX_distribution_path,treatment,"/",context_vec[context_n],"_CX_count.csv"))
      
      cntx_counts = cntx_counts[cntx_counts$type == ann_list[ann.l],]
      
      ### box plot
      cntrl_col = grep(paste0(control,"\\.[0-9]"), names(cntx_counts))
      trmnt_col = grep(paste0(treatment,"\\.[0-9]"), names(cntx_counts))
      
      ann_df = as.data.frame(t(cntx_counts[,c(cntrl_col,trmnt_col)]))
      names(ann_df)[1] = "ann_col"
      ann_df$tr = gsub("\\.[0-9]","", row.names(ann_df))
      ann_df$tr = factor(ann_df$tr, levels = c(control, treatment))

      
      
      proportion_value <- round(cntx_counts$ave.proportion, 3)
      max_y_value <- max(ann_df$ann_col, na.rm = TRUE)
      min_y_value <- min(ann_df$ann_col, na.rm = TRUE)
      
      ttest_p <- t.test(ann_df$ann_col ~ ann_df$tr)$p.value
      pValue_marks = ifelse(ttest_p < 0.001, "***", ifelse(ttest_p < 0.01, "**", ifelse(ttest_p < 0.05, "*", "")))
      
      # Create the ggplot boxplot
      p = ggplot(ann_df, aes(x = factor(tr), y = ann_col)) + 
        geom_boxplot(aes(fill = factor(tr)), ) + 
        scale_fill_manual(values = c("lightgray", color_vec[context_n])) +
        labs(
          title = paste0(context_vec[context_n]),  # Set the main title to the context value
          x = "", y = ""
        ) + 
        scale_y_continuous(#labels = scales::scientific,
                           limits = c(min_y_value, max_y_value+((max_y_value-min_y_value)*0.2)),
                           breaks = c(min_y_value,mean(c(min_y_value,max_y_value)),max_y_value)) +  # Format y-axis in scientific notation
        annotate("text", x = 1.5, y = max_y_value, 
                 label = paste0("X ", proportion_value), 
                 size = 4, hjust = 0.5, vjust = -1.75) +
        annotate("text", x = 1.5, y = max_y_value, 
                 label = pValue_marks,
                 size = 8, hjust = 0.5, vjust = 0.35) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 10),
          panel.background = element_blank(),    # Removes the background
          panel.border = element_rect(color = "black", fill = NA, size = 1), # Adds a box around the plot
          panel.grid.major = element_blank(),    # Removes major gridlines
          panel.grid.minor = element_blank(),    # Removes minor gridlines
          plot.background = element_blank(),     # Removes the plot background
          legend.position = "none",               # Removes the legend
          plot.margin = margin(t = 30, r = 5, b = 10, l = 0, unit = "pt") # r = ifelse(context_n == 15, 10, 0)
        )
      
      # Append the plot to the list
      plot_list[[context_n+1]] = p
    }
    
    ann_name_plot = ifelse(ann_list[ann.l] == "Retrotransposons",
                           "LINE/L1, LTR/Copya, LTR/Gypsy",
                           gsub("_"," ",ann_list[ann.l])
    )
    # save plot
    combined_plot <- plot_grid(plotlist = plot_list, ncol = 16,
                               labels = ann_name_plot,
                               label_size = 20, hjust = -0.025,
                               rel_widths = c(0.2,rep(1,15))) # first is yaxis_title, seconf is 14 contexts and the last is more becous is right border is farther
    print(combined_plot)
  }
  dev.off()
}
