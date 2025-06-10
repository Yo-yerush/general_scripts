library(dplyr)
library(ggplot2)

treatment = "mto1_vs_wt"

gff = rtracklayer::import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 gff3/TAIR10_GFF3_genes.gff")
chr_width = gff[which(gff$type == "chromosome")] %>% as.data.frame() %>% filter(width > 1e6) %>% .$width / 1e6


CG = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/",treatment,"/DMRs_CG_",treatment,".csv"))
CHG = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/",treatment,"/DMRs_CHG_",treatment,".csv"))
CHH = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/",treatment,"/DMRs_CHH_",treatment,".csv"))

###########################################################
#### count DMRs per 1Mbp
create_df <- function(x, n.chr, chr.w=chr_width) {
  chr = paste0("Chr",n.chr)
  x = x[x$seqnames == chr, "regionType"]
  x.up = length(x[x == "gain"]) / chr.w[n.chr]
  x.down = length(x[x == "loss"]) / chr.w[n.chr]
  
  return(c(chr,x.up,x.down))
}

chr_df_CG = data.frame(Chr = NA, Hyper = NA, Hypo = NA)
chr_df_CHG = data.frame(Chr = NA, Hyper = NA, Hypo = NA)
chr_df_CHH = data.frame(Chr = NA, Hyper = NA, Hypo = NA)

for (i in 1:5) {
  chr_df_CG[i,] = create_df(CG,i)
  chr_df_CHG[i,] = create_df(CHG,i)
  chr_df_CHH[i,] = create_df(CHH,i)
}

###########################################################
#### plot
plot_chr <- function(x.df,cntx) {
  df <- reshape2::melt(x.df, id.vars = "Chr", variable.name = "Type", value.name = "Count")
  df$Count = as.numeric(df$Count)
  
  ggplot(df, aes(x = Chr, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5, colour = "black") +
    scale_fill_manual(values = c("Hyper" = "#d96c6c", "Hypo" = "#6c96d9")) +
    theme_classic() +
    theme(legend.key.size = unit(0.4, "cm"), 
          text=element_text(family="serif"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, face="bold"),
          axis.text.y = element_text(face="bold", size=8, ),
          axis.title.y = element_text(size = 12, face="bold"),
          axis.line = element_blank(),
          axis.ticks = element_line(size = 0.75),
          axis.ticks.length.x = unit(-0.2, "cm"),
          axis.ticks.length.y = unit(-0.1, "cm"),
          legend.position = "none",
          title =  element_text(face="bold")
    ) +
    scale_y_continuous(breaks = seq(0, max(df$Count)*1.05, by = round(max(df$Count) / 5))) +
    labs(y = "DMRs Count per 1Mbp", x = "", fill = "Direction", title = cntx) +
    geom_rect(aes(xmin = 0.5, xmax = 5.5, ymin = 0, ymax = round(max(Count))*1.05),
              fill = "transparent", color = "black", size = 0.75)
}


legend_plot <- function(x.df) {
  df = data.frame(a = c(1,2), Type = c("Hyper","Hypo"), b = c(1,2))
  l.p = ggplot(df, aes(x = a, y = b, fill = Type)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("Hyper" = "#d96c6c", "Hypo" = "#6c96d9")) +
    theme_classic() +
    theme(#panel.spacing = unit(2, "lines"),
      legend.key.size = unit(0.4, "cm"), 
      text=element_text(family="serif"),
      legend.positio = "right") +
    labs(fill = "Direction")
  
  legend.p = suppressWarnings(cowplot::get_legend(l.p))
  cowplot::plot_grid(legend.p, nrow = 2)
}

CG_plot = plot_chr(chr_df_CG, "CG")
CHG_plot = plot_chr(chr_df_CHG, "CHG")
CHH_plot = plot_chr(chr_df_CHH, "CHH")

svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/",treatment,"_Chr_bar_plots.svg"),
    width = 6.25, height = 2.75)
cowplot::plot_grid(CG_plot, CHG_plot, CHH_plot, legend_plot(),
                   nrow = 1,
                   rel_widths = c(2, 2, 2, 1))
dev.off()

