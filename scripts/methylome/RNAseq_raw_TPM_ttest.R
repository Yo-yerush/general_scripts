library(plyr)
library(dplyr)
library(corrplot)
library(ggpubr)
library(ggplot2)
library(lmerTest)
library(tidyr)

for (treatment in c("mto1","mto3","35s")) {
  
  path_2_save = "P:/yonatan/methionine/rnaseq_23/raw_TPM_ttest/"
  
  if(treatment == "mto1") {
    mto_rnaseq_names = c("met14","met15","met16")
    #wt_rnaseq_names = c("met20","met22")
    wt_rnaseq_names = c("met20","met21","met22")
    
  } else  if (treatment == "mto3") {
    mto_rnaseq_names = c("met17","met18","met19")
    #wt_rnaseq_names = c("met20","met22")
    wt_rnaseq_names = c("met20","met21","met22")
    
  } else  if (treatment == "35s") {
    mto_rnaseq_names = c("met24","met25","met26")
    #wt_rnaseq_names = c("met20","met22")
    wt_rnaseq_names = c("met27","met28","met29")
    
  }
  
  
  #Load the RNAseq '*.gene.results' files (output from RSEM pipeline)
  read_and_select <- function(filename) {
    filepath <- paste0("P:/yonatan/methionine/rnaseq_23/gene.results.files/", filename, ".genes.results")
    df <- read.table(filepath, header = TRUE)
    df <- df %>% select(gene_id, TPM)
    names(df)[2] = filename
    return(df)
  }
  RNA.0 <- suppressMessages(bind_cols(lapply(c(wt_rnaseq_names,mto_rnaseq_names), read_and_select)))
  names(RNA.0)[1] = "transcript_id"
  RNA.0 = RNA.0[,-grep("gene.id",names(RNA.0))]
  

  #### Gene expression ####
  refseq_2_tair = read.csv("P:/TAIR10.1/TAIR_2_refseq_transcript.csv")[,1:2]
  RNA = merge.data.frame(refseq_2_tair, RNA.0, by = "transcript_id", all.y = T)
  # edit table and names
  mto_names_pos = grep(paste(mto_rnaseq_names, collapse = "|"), names(RNA))
  names(RNA)[mto_names_pos] = paste0(treatment,".",1:3,"_TPM")
  
  wt_names_pos = grep(paste(wt_rnaseq_names, collapse = "|"), names(RNA))
  names(RNA)[wt_names_pos] = paste0("wt.",1:3,"_TPM")
  #names(RNA)[wt_names_pos] = paste0("wt.",1:2,"_TPM")
  
  RNA = RNA[,c(1,2,wt_names_pos,mto_names_pos)]
  
  for (i in 1:nrow(RNA)) {
    mto_l_values = as.numeric(RNA[i, mto_names_pos])
    wt_l_values = as.numeric(RNA[i, wt_names_pos])
    
    RNA$log2FC[i] = log2(mean(mto_l_values) / mean(wt_l_values))
    try({RNA$pValue[i] = t.test(mto_l_values, wt_l_values, var.equal = T)$p.value})
  }
  
  RNA = RNA[order(RNA$pValue),]
  
  write.csv(RNA, paste0(path_2_save,treatment,"_TPM_ttest.csv"), row.names = F)
}
