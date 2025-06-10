source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/revigo_function.R")

path2results = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/revigo_results/"
methylome_GO_path = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/"
RNAseq_GO_path = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/GO_analysis/results/"

#######
# darwin
source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/revigo_function.R")

path2results = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/revigo_results/"
methylome_GO_path = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/"
RNAseq_GO_path = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/GO_analysis/results/"
#######

##############################
# save scatter plot function
scatter_plot_run_function <- function(x,x.file, x.file.suffix) {
  # save scatter plot
  svg(paste0(x.file,"up_ScatterPlot",x.file.suffix,".svg"), width = 5, height = 2.38, family = "serif")
  print(x$scatterPlot_up)
  dev.off()
  
  svg(paste0(x.file,"down_ScatterPlot",x.file.suffix,".svg"), width = 5, height = 2.38, family = "serif")
  print(x$scatterPlot_down)
  dev.off()
  
  # save scatter plot legend
  svg(paste0(x.file,"ScatterPlot_legend",x.file.suffix,".svg"), width = 4.32, height = 10, family = "serif")
  par(mar = c(0.1, 0.1, 1.1, 0.1))
  plot.new()
  legend("topleft",
         legend = x$scatterPlot_legend$parentTerm,
         pch = 21,
         pt.cex = 1.5,
         pt.bg = paste0(x$scatterPlot_legend$hex_col,75),
         col = "black", 
         cex = 0.75, 
         bty = "n",
         x.intersp = 0.8,
         y.intersp = 1)
  title(main = "Parent Term", cex.main = 1, adj = 0.01)
  dev.off()
  
  svg(paste0(x.file,"up_treeMap",x.file.suffix,".svg"), width = 5, height = 5, family = "serif")
  treemapPlot(x$treemapPlot_up)
  dev.off()
  
  svg(paste0(x.file,"down_treeMap",x.file.suffix,".svg"), width = 5, height = 5, family = "serif")
  treemapPlot(x$treemapPlot_down)
  dev.off()
  
  # save scatter plot data frame
  write.csv(x$revigo_df, paste0(x.file,"ScatterData",x.file.suffix,".csv"), row.names = F) 
}

##############################

############### RNAseq results
for (treat.l in c("mto1_vs_wt","mto3_vs_wt","dCGS_vs_EV", "SSE_H_vs_EV", "SSE_L_vs_EV")) {
  
  # manage final directory
  dir.create(paste0(path2results,treat.l), showWarnings = F)
  path2results_rna = paste0(path2results,treat.l,"/RNAseq/")
  dir.create(path2results_rna, showWarnings = F)
  
  for (type.l in c("BP","MF","CC")) {
    
    path2results_rna_type = paste0(path2results_rna,type.l,"/")
    dir.create(path2results_rna_type)
    tryCatch({
      GO_file_up = read.csv(paste0(RNAseq_GO_path,treat.l,"/",type.l,"/topGO_",type.l,"_up_",treat.l,".csv"))
      GO_file_down = read.csv(paste0(RNAseq_GO_path,treat.l,"/",type.l,"/topGO_",type.l,"_down_",treat.l,".csv"))
      
      #####################################
      
      scatterPlot_fun = Revigo_plots(GO_df_up = GO_file_up,
                                     GO_df_down = GO_file_down,
                                     #treatment = treat.l,
                                     GO_type = type.l)
      file_suffix = paste0("_",type.l,"_",treat.l,"_RNAseq")
      
      # save scatter plot
      scatter_plot_run_function(scatterPlot_fun, path2results_rna_type, file_suffix)
      
    }, error = function(e){message(paste0("\n**\tfail: ",treat.l,"_",type.l,"\t**\n"))
    })
  }
}

##############################
############### Methylome results
for (treat.l in c("mto1_vs_wt","mto3_vs_wt","dCGS_vs_EV", "SSE_H_vs_EV", "SSE_L_vs_EV")) {
  # manage final directory
  dir.create(paste0(path2results,treat.l), showWarnings = F)
  path2results_meth = paste0(path2results,treat.l,"/Methylome/")
  dir.create(path2results_meth, showWarnings = F)
  
  for (ann.l in c("Genes","Promoters", "CDS", "Introns")) {
    path2results_meth_2 = paste0(path2results_meth,ann.l,"/")
    dir.create(path2results_meth_2, showWarnings = F)
    
    for (context.l in c("CG","CHG","CHH","all")) {
      path2results_meth_3 = paste0(path2results_meth_2,context.l,"/")
      dir.create(path2results_meth_3, showWarnings = F)
      
      for (type.l in c("BP","MF","CC")) {
        path2results_meth_type = paste0(path2results_meth_3,type.l,"/")
        dir.create(path2results_meth_type)
        
        tryCatch({
          if (context.l != "all") {
            GO_file_up = read.csv(paste0(methylome_GO_path,treat.l,"/GO_analysis/",context.l,"/",ann.l,"/",type.l,".",ann.l,".",context.l,".gain.",treat.l,".topGO.csv"))
            GO_file_down = read.csv(paste0(methylome_GO_path,treat.l,"/GO_analysis/",context.l,"/",ann.l,"/",type.l,".",ann.l,".",context.l,".loss.",treat.l,".topGO.csv"))
          } else {
            GO_file_up = rbind(read.csv(paste0(methylome_GO_path,treat.l,"/GO_analysis/CG/",ann.l,"/",type.l,".",ann.l,".CG.gain.",treat.l,".topGO.csv")),
                               read.csv(paste0(methylome_GO_path,treat.l,"/GO_analysis/CHG/",ann.l,"/",type.l,".",ann.l,".CHG.gain.",treat.l,".topGO.csv")),
                               read.csv(paste0(methylome_GO_path,treat.l,"/GO_analysis/CHH/",ann.l,"/",type.l,".",ann.l,".CHH.gain.",treat.l,".topGO.csv")))
            GO_file_down = rbind(read.csv(paste0(methylome_GO_path,treat.l,"/GO_analysis/CG/",ann.l,"/",type.l,".",ann.l,".CG.loss.",treat.l,".topGO.csv")),
                                 read.csv(paste0(methylome_GO_path,treat.l,"/GO_analysis/CHG/",ann.l,"/",type.l,".",ann.l,".CHG.loss.",treat.l,".topGO.csv")),
                                 read.csv(paste0(methylome_GO_path,treat.l,"/GO_analysis/CHH/",ann.l,"/",type.l,".",ann.l,".CHH.loss.",treat.l,".topGO.csv")))
          }
          names(GO_file_up) = gsub("Fisher", "pValue", names(GO_file_up))
          names(GO_file_down) = gsub("Fisher", "pValue", names(GO_file_down))
          
          #####################################
          
          scatterPlot_fun = Revigo_plots(GO_df_up = GO_file_up,
                                         GO_df_down = GO_file_down,
                                         #treatment = treat.l,
                                         GO_type = type.l)
          file_suffix = paste0("_",type.l,"_",treat.l,"_Methylome")
          
          # save scatter plot
          scatter_plot_run_function(scatterPlot_fun, path2results_meth_type, file_suffix)
          
        }, error = function(e){message(paste0("\n**\tfail: ",treat.l,"_",ann.l,"_",context.l,"_",type.l,"\t**\n"))
        })
      }
    }
  }
}