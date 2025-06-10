if (F) {
  experiment = "test"
  name1 = "mto1_vs_wt"
  name2 = "mto3_vs_wt"
  name3 = mto1_rnaseq
  name4 = mto3_rnaseq
  name5 = "sse_high_vs_sse_ev"
  DMRcaller.dir = "P:/yonatan/methionine/BSseq_pipeline/DMRcaller/DMRcaller_res/"
  context = "CG"
  save.files = NULL
}

venn.di = function(experiment,
                   name1, # "A_vs_B"
                   name2, # "C_vs_D"
                   name3 = NA,
                   name4 = NA,
                   name5 = NA,
                   DMRcaller.dir = "P:/yonatan/methionine/BSseq_results_100823/", # path of DMRcaller resulrs 
                   save.files = "P:/yonatan/methionine/BSseq_results_100823/DMRs_venn/")
{
  
  library(dplyr)
  library(ggplot2)
  library(VennDiagram)
  library("ggVennDiagram")
  library(RColorBrewer)
  
  ###
  for (context in c("CG","CHG","CHH")) {
    
    dir.create(paste0(save.files,experiment))
    
    if (file.exists(paste0(save.files,experiment,"/",experiment,"_",context,"_merged.csv")) == T) {
      stop("change <experiment> name")
    }
    
    
    # load the data
    #column_keep = c("seqnames","start","end","cytosinesCount","context","pValue","regionType")
    column_keep = c("seqnames","start","end","regionType")
    
    set1 = read.csv(paste0(DMRcaller.dir,name1,"/DMRs_",context,"_",name1,".csv"))[,column_keep]
    set2 = read.csv(paste0(DMRcaller.dir,name2,"/DMRs_",context,"_",name2,".csv"))[,column_keep]
    
    # edit location column
    set1$DMR_loc = paste0(set1[,1],":",set1[,2],"-",set1[,3])
    set2$DMR_loc = paste0(set2[,1],":",set2[,2],"-",set2[,3])
    set1 = set1 %>% select(DMR_loc, everything()) %>% select(-c(2:4))
    set2 = set2 %>% select(DMR_loc, everything()) %>% select(-c(2:4))
    
    # load for 3, 4 and 5 sets
    if(is.na(name5) == F) {
      set3 = read.csv(paste0(DMRcaller.dir,name3,"/DMRs_",context,"_",name3,".csv"))[,column_keep]
      set4 = read.csv(paste0(DMRcaller.dir,name4,"/DMRs_",context,"_",name4,".csv"))[,column_keep]
      set5 = read.csv(paste0(DMRcaller.dir,name5,"/DMRs_",context,"_",name5,".csv"))[,column_keep]
      #
      set3$DMR_loc = paste0(set3[,1],":",set3[,2],"-",set3[,3])
      set4$DMR_loc = paste0(set4[,1],":",set4[,2],"-",set4[,3])
      set5$DMR_loc = paste0(set5[,1],":",set5[,2],"-",set5[,3])
      set3 = set3 %>% select(DMR_loc, everything()) %>% select(-c(2:4))
      set4 = set4 %>% select(DMR_loc, everything()) %>% select(-c(2:4))
      set5 = set5 %>% select(DMR_loc, everything()) %>% select(-c(2:4))
      #
      category.names = c(name1, name2, name3, name4, name5)
    } else if(is.na(name4) == F) {
      set3 = read.csv(paste0(DMRcaller.dir,name3,"/DMRs_",context,"_",name3,".csv"))[,column_keep]
      set4 = read.csv(paste0(DMRcaller.dir,name4,"/DMRs_",context,"_",name4,".csv"))[,column_keep]
      #
      set3$DMR_loc = paste0(set3[,1],":",set3[,2],"-",set3[,3])
      set4$DMR_loc = paste0(set4[,1],":",set4[,2],"-",set4[,3])
      set3 = set3 %>% select(DMR_loc, everything()) %>% select(-c(2:4))
      set4 = set4 %>% select(DMR_loc, everything()) %>% select(-c(2:4))
      #
      category.names = c(name1, name2, name3, name4)
    } else if (is.na(name3) == F) {
      set3 = read.csv(paste0(DMRcaller.dir,name3,"/DMRs_",context,"_",name3,".csv"))[,column_keep]
      #
      set3$DMR_loc = paste0(set3[,1],":",set3[,2],"-",set3[,3])
      set3 = set3 %>% select(DMR_loc, everything()) %>% select(-c(2:4))
      #
      category.names = c(name1, name2, name3)
    } else {
      category.names = c(name1, name2)
    }
    
    
    if(is.na(name5) == F) {
      x = list(set1$DMR_loc, set2$DMR_loc, set3$DMR_loc, set4$DMR_loc, set5$DMR_loc)
    } else if(is.na(name4) == F) {
      x = list(set1$DMR_loc, set2$DMR_loc, set3$DMR_loc, set4$DMR_loc)
    } else if (is.na(name3) == F) {
      x = list(set1$DMR_loc, set2$DMR_loc, set3$DMR_loc)
    } else {
      x = list(set1$DMR_loc, set2$DMR_loc)
    }
    
    ###
    
    # venn diagram position and resolution
    if(is.na(name5) == F) {
      category.position = c(0,0,234,160,0)
      resolution = 250
      cex = 0.3
    } else if(is.na(name4) == F) {
      category.position = c(-10, 5, 0, 0)
      resolution = 250
      cex = 0.5
    } else if (is.na(name3) == F) {
      category.position = c(-27, 27, 180)
      resolution = 300
      cex = 0.5
    } else {
      category.position = c(0, 0)
      resolution = 300
      cex = 0.5
    }
    
    #fill = c("#440154ff", '#21908dff', '#fde725ff')
    venn_colors = c("red", "blue","yellow", "green", "purple")
    
    venn.diagram(
      x = x,
      category.names = gsub("_"," ",category.names),
      filename = paste0(save.files,experiment,"/",experiment,"_",context,"_VennDiagram.png"),
      disable.logging = T,
      output = TRUE ,
      imagetype="png" ,
      height = 480 , 
      width = 480 , 
      resolution = resolution,
      compression = "lzw",
      lwd = 1,
      fill = venn_colors[1:length(x)],
      alpha = rep(0.3, length(x)),
      cex = cex,
      fontfamily = "sans",
      cat.cex = 0.4,
      cat.default.pos = "outer",
      cat.pos = category.position,
      cat.fontface = 2,
      cat.fontfamily = "sans"
      #    cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
      #    col=venn_colors,
      #    rotation = 1
    )
    
    
    ###
    
    sets_merge = merge.data.frame(set1,set2, by = "DMR_loc")
    
    if(is.na(name5) == F) {
      sets_merge.2 = merge.data.frame(sets_merge, set3, by = "DMR_loc")
      sets_merge.3 = merge.data.frame(sets_merge.2, set4, by = "DMR_loc")
      merge.file = merge.data.frame(sets_merge.3, set5, by = "DMR_loc")
      colnames(merge.file)[2:6] = category.names
      
    } else if(is.na(name4) == F) {
      sets_merge.2 = merge.data.frame(sets_merge, set3, by = "DMR_loc")
      merge.file = merge.data.frame(sets_merge.2, set4, by = "DMR_loc")
      colnames(merge.file)[2:5] = category.names
      
    } else if (is.na(name3) == F) {
      merge.file = merge.data.frame(sets_merge, set3, by = "DMR_loc")
      colnames(merge.file)[2:4] = category.names
      
    } else {
      merge.file = sets_merge
      colnames(merge.file)[2:3] = category.names
    }
    
    if (is.null(save.files) == F) {
      write.csv(merge.file, paste0(save.files,experiment,"/",experiment,"_",context,"_merged.csv"), row.names=FALSE)
    } else {merge.file}
    
  }
}