source("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/scripts/revigo_function.R")
path2results = "C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/revigo_results/"

##############################
############### RNAseq results
for (treat.l in c("mto1_vs_wt")) {
  
  # manage final directory
  dir.create(paste0(path2results,treat.l), showWarnings = F)
  path2results_rna = paste0(path2results,treat.l,"/RNAseq/")
  dir.create(path2results_rna, showWarnings = F)
  
  for (type.l in c("BP","MF","CC")) {
    for (direction.l in c("up","down")) {
      
      tryCatch({
        GO_file_up = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/GO_analysis/results/",treat.l,"/",type.l,"/topGO_",type.l,"_up_",treat.l,".csv"))
        GO_file_down = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/GO_analysis/results/",treat.l,"/",type.l,"/topGO_",type.l,"_down_",treat.l,".csv"))
        
        scatterPlot_fun = Revigo_plots(GO_df_up = GO_file_up,
                                       GO_df_down = GO_file_down,
                                       treatment = treat.l,
                                       GO_type = type.l)
        
        file_name = paste0(path2results_rna,treat.l,"_",type.l,"_Revigo_")
        
        # save scatter plot
        svg(paste0(file_name,"ScatterPlot.svg"), width = 5.04, height = 4.17, family = "serif")
        print(scatterPlot_fun$scatterPlot)
        dev.off()
        
        # save scatter plot data frame
        write.csv(scatterPlot_fun$revigo_df, paste0(file_name,"ScatterData.csv"), row.names = F)  
        
        # save tree map
        svg(paste0(file_name,"treeMap.svg"), width = 5, height = 5, family = "serif")
        treemapPlot(scatterPlot_fun$treemapPlot_vec)
        dev.off()
        
      }, error = function(e){message(paste0("\n**\tfail: ",treat.l,"_",type.l,"\t**\n"))
      })
    }
  }
}

##############################
############### Methylome results
for (treat.l in c("mto1_vs_wt","mto3_vs_wt","dCGS_vs_EV")) {
  # manage final directory
  dir.create(paste0(path2results,treat.l), showWarnings = F)
  path2results_meth = paste0(path2results,treat.l,"/Methylome/")
  dir.create(path2results_meth, showWarnings = F)
  
  for (ann.l in c("Genes","Promoters","Introns","TEG")) {
    path2results_meth_2 = paste0(path2results_meth,ann.l,"/")
    dir.create(path2results_meth_2, showWarnings = F)
    
    for (context.l in c("CG","CHG","CHH","all")) {
      path2results_meth_3 = paste0(path2results_meth_2,context.l,"/")
      dir.create(path2results_meth_3, showWarnings = F)
      
      for (type.l in c("BP","MF","CC")) {
        for (direction.l in c("gain","loss")) {
          
          tryCatch({
            if (context.l != "all") {
              GO_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_191223/",treat.l,"/GO_analysis/",context.l,"/",ann.l,"/",type.l,".",ann.l,".",context.l,".",direction.l,".",treat.l,".topGO.csv"))
            } else {
              GO_file = rbind(read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_191223/",treat.l,"/GO_analysis/CG/",ann.l,"/",type.l,".",ann.l,".CG.",direction.l,".",treat.l,".topGO.csv")),
                              read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_191223/",treat.l,"/GO_analysis/CHG/",ann.l,"/",type.l,".",ann.l,".CHG.",direction.l,".",treat.l,".topGO.csv")),
                              read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_191223/",treat.l,"/GO_analysis/CHH/",ann.l,"/",type.l,".",ann.l,".CHH.",direction.l,".",treat.l,".topGO.csv")))
            }
            names(GO_file) = gsub("Fisher", "pValue", names(GO_file))
            
            scatterPlot_fun = Revigo_plots(GO_df = GO_file,
                                           treatment = treat.l,
                                           direction = direction.l,
                                           GO_type = type.l)
            
            file_name = paste0(path2results_meth_3,treat.l,"_",ann.l,"_",context.l,"_",type.l,"_Revigo_")
            
            # save scatter plot
            svg(paste0(file_name,"ScatterPlot.svg"), width = 5.04, height = 4.17, family = "serif")
            print(scatterPlot_fun$scatterPlot)
            dev.off()
            
            # save scatter plot data frame
            write.csv(scatterPlot_fun$revigo_df, paste0(file_name,"ScatterData.csv"), row.names = F) 
            
            # save tree map
            svg(paste0(file_name,"treeMap.svg"), width = 5, height = 5, family = "serif")
            treemapPlot(scatterPlot_fun$treemapPlot_vec)
            dev.off()
            
          }, error = function(e){message(paste0("\n**\tfail: ",treat.l,"_",ann.l,"_",context.l,"_",type.l,"\t**\n"))
          })
        }
      }
    }
  }
}

