venn.di = function(experiment, up.down, name1 = "A_VS_B", name2 = "B_VS_C", name3 = NA, name4 = NA, 
                   merge.file.number, experiment.2 = NA, up.down.2 = up.down, path = getwd()) {

  library(ggplot2)
  library(VennDiagram)
  library("ggVennDiagram")
  library(RColorBrewer)
  
  if (is.na(experiment.2) == F) {
    path.2 = paste(path,experiment.2,"statistics/", sep = "/")
  }
  path = paste(path,experiment,"statistics/", sep = "/")
  
  ###
  if (file.exists(paste(experiment,merge.file.number,up.down,"merge", "csv", sep = ".")) == T) {
    stop("change <merge.file.number>")
  }
  
  set1 = read.csv(paste(paste(path,name1,"/",name1, sep = ""),"FC",up.down,"DE.csv",sep = "."))
  
  if (is.na(experiment.2) == F) {
    set2 = read.csv(paste(paste(path.2,name2,"/",name2, sep = ""),"FC",up.down.2,"DE.csv",sep = "."))
  } else {
    set2 = read.csv(paste(paste(path,name2,"/",name2, sep = ""),"FC",up.down,"DE.csv",sep = "."))
  }
  
  if(is.na(name4) == F) {
    category.names = c(paste(name1,"log2FC",sep = "_"), paste(name2,"log2FC",sep = "_"), paste(name3,"log2FC",sep = "_"), paste(name4,"log2FC",sep = "_"))
  } else if (is.na(name3) == F) {
    category.names = c(paste(name1,"log2FC",sep = "_"), paste(name2,"log2FC",sep = "_"), paste(name3,"log2FC",sep = "_"))
  } else {
    category.names = c(paste(name1,"log2FC",sep = "_"), paste(name2,"log2FC",sep = "_"))
  }
  
  
  if(is.na(name4) == F) {
    set3 = read.csv(paste(paste(path,name3,"/",name3, sep = ""),"FC",up.down,"DE.csv",sep = "."))
    set4 = read.csv(paste(paste(path,name4,"/",name4, sep = ""),"FC",up.down,"DE.csv",sep = "."))
  } else if (is.na(name3) == F) {
    set3 = read.csv(paste(paste(path,name3,"/",name3, sep = ""),"FC",up.down,"DE.csv",sep = "."))
  } else {
    NA
  }
  
  
  if(is.na(name4) == F) {
    x = list(set1$transcript_id, set2$transcript_id, set3$transcript_id, set4$transcript_id)
  } else if (is.na(name3) == F) {
    x = list(set1$transcript_id, set2$transcript_id, set3$transcript_id)
  } else {
    x = list(set1$transcript_id, set2$transcript_id)
  }
  
  ###
  
  setwd(path)
  
  plot.venn = ggVennDiagram(x, label_alpha = 0, category.names = category.names) + 
    scale_fill_gradient(low="white",high = "red")
  
  if (is.na(experiment.2) == F) {
    ggsave(paste(experiment,up.down,experiment.2,up.down.2,merge.file.number,"venn.diagram","png", sep = "."), plot = plot.venn)
  } else {
    ggsave(paste(experiment,up.down,merge.file.number,"venn.diagram","png", sep = "."), plot = plot.venn)
  }
  
  ###
  
  setwd(path)
  
  sets_merge = merge.data.frame(set1,set2, by = "transcript_id")
  
  if(is.na(name4) == F) {
    sets_merge.2 = merge.data.frame(sets_merge, set3, by = "transcript_id")
    sets_merge.3 = merge.data.frame(sets_merge.2, set4, by = "transcript_id")
  } else if (is.na(name3) == F) {
    sets_merge.2 = merge.data.frame(sets_merge, set3, by = "transcript_id")
  } else {
    NA
  }
  
  
  if(is.na(name4) == F) {
    merge.file = sets_merge.3[,-c(5:9,12:16,19:23)]
  } else if (is.na(name3) == F) {
    merge.file = sets_merge.2[,-c(5:9,12:16)]
  } else {
    merge.file = sets_merge[,-c(5:9)]
  }
  
  
  if(is.na(name4) == F) {
    merge.file.col.names = c(3,5,7,9)
  } else if (is.na(name3) == F) {
    merge.file.col.names = c(3,5,7)
  } else {
    merge.file.col.names = c(3,5)
  }
  
  colnames(merge.file)[merge.file.col.names] = category.names
  
  if (is.na(experiment.2) == F) {
    write.csv(merge.file, paste(experiment,up.down,experiment.2,up.down.2,merge.file.number,"merge", "csv", sep = "."), row.names=FALSE)
  } else {
    write.csv(merge.file, paste(experiment,up.down,merge.file.number,"merge", "csv", sep = "."), row.names=FALSE)
  }
}