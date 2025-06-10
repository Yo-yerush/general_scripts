library(DMRcaller)
library(dplyr)
library(parallel)

samples_mto1 = read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_mto1.txt", sep = "\t", header = F)
samples_mto3 = read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_mto3.txt", sep = "\t", header = F)
samples_35s = read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_35s.txt", sep = "\t", header = F)

wt_path = samples_mto1[1:2,2]
mto1_path = samples_mto1[3:5,2]
mto3_path = samples_mto3[3:5,2]
ev_path = samples_35s[1:3,2]
p35s_path = samples_35s[4:6,2]

source("/home/yoyerush/yo/methylome_pipeline/Methylome.At/Methylome.At_scripts/load_replicates.R")
load_wt = load_replicates(wt_path, n.cores = 3)
load_mto1 = load_replicates(mto1_path, n.cores = 3)
load_mto3 = load_replicates(mto3_path, n.cores = 3)
load_ev = load_replicates(ev_path, n.cores = 3)
load_35s = load_replicates(p35s_path, n.cores = 3)

# trimm seqs objects
trimm_Chr <- function(gr_obj) {
  remove_seqnames = c("NC_000932.1","NC_037304.1")
  gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
  seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
  return(sort(gr_obj))
}
meth_wt = trimm_Chr(load_wt$methylationDataReplicates)
meth_mto1 = trimm_Chr(load_mto1$methylationDataReplicates)
meth_mto3 = trimm_Chr(load_mto3$methylationDataReplicates)
meth_ev = trimm_Chr(load_ev$methylationDataReplicates)
meth_35s = trimm_Chr(load_35s$methylationDataReplicates)

total_meth_df <- function(meth_ratio_var1) {
  
  meth_ratio_var1$readsT1 = (meth_ratio_var1$readsM1/meth_ratio_var1$readsN1)*100
  meth_ratio_var1$readsT2 = (meth_ratio_var1$readsM2/meth_ratio_var1$readsN2)*100
  if (length(grep("readsM3",names(meth_ratio_var1@elementMetadata))) == 1) {
    meth_ratio_var1$readsT3 = (meth_ratio_var1$readsM3/meth_ratio_var1$readsN3)*100
    var1_t3 = T
  } else {var1_t3 = F}
  
  
  ############ average for all ration separately 
  var1_total_CG = meth_ratio_var1[which(meth_ratio_var1$context == "CG")]
  var1_total_CHG = meth_ratio_var1[which(meth_ratio_var1$context == "CHG")]
  var1_total_CHH = meth_ratio_var1[which(meth_ratio_var1$context == "CHH")]
  if (var1_t3) {
    var1_CG = c(mean(var1_total_CG@elementMetadata$readsT1, na.rm = T),
                mean(var1_total_CG@elementMetadata$readsT2, na.rm = T),
                mean(var1_total_CG@elementMetadata$readsT3, na.rm = T))
    var1_CHG = c(mean(var1_total_CHG@elementMetadata$readsT1, na.rm = T),
                 mean(var1_total_CHG@elementMetadata$readsT2, na.rm = T),
                 mean(var1_total_CHG@elementMetadata$readsT3, na.rm = T))
    var1_CHH = c(mean(var1_total_CHH@elementMetadata$readsT1, na.rm = T),
                 mean(var1_total_CHH@elementMetadata$readsT2, na.rm = T),
                 mean(var1_total_CHH@elementMetadata$readsT3, na.rm = T))
  } else {
    var1_CG = c(mean(var1_total_CG@elementMetadata$readsT1, na.rm = T),
                mean(var1_total_CG@elementMetadata$readsT2, na.rm = T))
    var1_CHG = c(mean(var1_total_CHG@elementMetadata$readsT1, na.rm = T),
                 mean(var1_total_CHG@elementMetadata$readsT2, na.rm = T))
    var1_CHH = c(mean(var1_total_CHH@elementMetadata$readsT1, na.rm = T),
                 mean(var1_total_CHH@elementMetadata$readsT2, na.rm = T))
  }
  
  
  meth_plot_df = data.frame(type = c("CG", "CHG", "CHH"),
                            levels = c(mean(var1_CG),mean(var1_CHG),mean(var1_CHH)),
                            SD = c(sd(var1_CG),sd(var1_CHG),sd(var1_CHH)))
  return(meth_plot_df)   
}

xl_list = list(
  wt = total_meth_df(meth_wt),
  mto1 = total_meth_df(meth_mto1),
  mto3 = total_meth_df(meth_mto3),
  ev_35s = total_meth_df(meth_ev),
  p35s = total_meth_df(meth_35s))

writexl::write_xlsx(xl_list, "/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/total_meth_df.xlsx")