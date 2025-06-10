library(GenomicRanges)
library(DMRcaller)
library(parallel)
library(ggplot2)
library(dplyr)
library(ggbreak)

treatment = "mto1"

sample_df = read.table(paste0("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_",treatment,".txt"))

TE_df = read.csv("/home/yoyerush/yo/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt", sep = "\t")

TE_df$Transposon_ID = TE_df$Transposon_Name
TE_df$Transposon_Name = gsub(paste0("AT1TE.*"), "NC_003070.9",TE_df$Transposon_Name)
TE_df$Transposon_Name = gsub(paste0("AT2TE.*"), "NC_003071.7",TE_df$Transposon_Name)
TE_df$Transposon_Name = gsub(paste0("AT3TE.*"), "NC_003074.8",TE_df$Transposon_Name)
TE_df$Transposon_Name = gsub(paste0("AT4TE.*"), "NC_003075.7",TE_df$Transposon_Name)
TE_df$Transposon_Name = gsub(paste0("AT5TE.*"), "NC_003076.8",TE_df$Transposon_Name)

TE_df = TE_df[,-c(5:6)]
names(TE_df)[1:4] = c("seqnames","strand","start","end")
TE_df$strand = ifelse(TE_df$strand == "true","+","-")
TE_df$strand = "*"
TE_gr = makeGRangesFromDataFrame(TE_df, keep.extra.columns = T) # dont need families characterization here

TE_gr$size <- width(TE_gr)
long_TEs = TE_gr[which(TE_gr$size > 4000)]
short_TEs = TE_gr[which(TE_gr$size < 500)]

##########################################

read_function <- function(file_path) {
  readBismark(file_path)
}

# Use mclapply to read files in parallel
read_CX <- mclapply(sample_df[,2], read_function, mc.cores = nrow(sample_df))
gr_list <- GRangesList(read_CX)
names(gr_list) = ave(sample_df[,1], sample_df[,1], FUN = function(x) paste0(x, ".", seq_along(x)))

##########################################

# create data frame of treatments and contexts
TE_ave_df <- data.frame(treatment = names(gr_list), CG = NA, CHG = NA, CHH = NA) 
TE_cor_df <- data.frame(treatment = names(gr_list), CG = NA, CHG = NA, CHH = NA)  

long_ave_df <- data.frame(treatment = names(gr_list), CG = NA, CHG = NA, CHH = NA) 
long_cor_df <- data.frame(treatment = names(gr_list), CG = NA, CHG = NA, CHH = NA)  

short_ave_df <- data.frame(treatment = names(gr_list), CG = NA, CHG = NA, CHH = NA) 
short_cor_df <- data.frame(treatment = names(gr_list), CG = NA, CHG = NA, CHH = NA)  

for (context in c("CG","CHG","CHH")) {
  
  ### count CX from annotations
  overlap_TE2trnt <- function(meth_TE, meth_cx) {
    
    meth_cx = meth_cx[which(meth_cx$context == context)]
    
    meth_cx = meth_cx[which(meth_cx$readsN >= 6)]
    meth_cx$proportion = meth_cx$readsM / meth_cx$readsN
    #meth_cx = meth_cx[which(meth_cx$proportion >= 0.4)]
    
    # find overlaps for all CX with TE annotation file
    ## var1 (control)
    m1 <- findOverlaps(meth_cx, meth_TE)
    CX_trnt <- meth_cx[queryHits(m1)]
    TE_trnt <- meth_TE[subjectHits(m1)]
    mcols(CX_trnt) <- cbind(mcols(CX_trnt), mcols(TE_trnt))
    
    mean_proportion = mean(CX_trnt$proportion)
    spearman_cor = cor(CX_trnt$proportion, CX_trnt$size, method = "spearman")
    
    return(list(mean_proportion = mean_proportion,
                spearman_cor = spearman_cor))
  }
  
  
  # run function on average and correlation data frames
  for (trnt.n in 1:length(names(gr_list))) {
    
    TE_ave_df[trnt.n,context] = overlap_TE2trnt(TE_gr, gr_list[[trnt.n]])$mean_proportion
    TE_cor_df[trnt.n,context] = overlap_TE2trnt(TE_gr, gr_list[[trnt.n]])$spearman_cor
    
    long_ave_df[trnt.n,context] = overlap_TE2trnt(long_TEs, gr_list[[trnt.n]])$mean_proportion
    long_cor_df[trnt.n,context] = overlap_TE2trnt(long_TEs, gr_list[[trnt.n]])$spearman_cor
    
    short_ave_df[trnt.n,context] = overlap_TE2trnt(short_TEs, gr_list[[trnt.n]])$mean_proportion
    short_cor_df[trnt.n,context] = overlap_TE2trnt(short_TEs, gr_list[[trnt.n]])$spearman_cor
    
    cat(">")
  }
  cat("\n")
}
##################################################################


################
## put this shit inside the loop
################

# t.test function for last row
ttest_fun <- function(x) {
  print(c("pValue",
          round(t.test(x$CG[1:2], x$CG[3:5])$p.value, 4),
          round(t.test(x$CHG[1:2], x$CHG[3:5])$p.value, 4),
          round(t.test(x$CHH[1:2], x$CHH[3:5])$p.value, 4)))
  
}

TE_ave_df[6,] = ttest_fun(TE_ave_df)
long_ave_df[6,] = ttest_fun(long_ave_df)
short_ave_df[6,] = ttest_fun(short_ave_df)

TE_cor_df[6,] = ttest_fun(TE_cor_df)
long_cor_df[6,] = ttest_fun(long_cor_df)
short_cor_df[6,] = ttest_fun(short_cor_df)

write.csv(TE_ave_df, "/home/yoyerush/yo/methylome_pipeline/TEs_methylation_ave_cor_LongShort/all_TEs_ave_df.csv", row.names = F)
write.csv(long_ave_df, "/home/yoyerush/yo/methylome_pipeline/TEs_methylation_ave_cor_LongShort/long_ave_df.csv", row.names = F)
write.csv(short_ave_df, "/home/yoyerush/yo/methylome_pipeline/TEs_methylation_ave_cor_LongShort/short_ave_df.csv", row.names = F)
write.csv(TE_cor_df, "/home/yoyerush/yo/methylome_pipeline/TEs_methylation_ave_cor_LongShort/all_TEs_cor_df.csv", row.names = F)
write.csv(long_cor_df, "/home/yoyerush/yo/methylome_pipeline/TEs_methylation_ave_cor_LongShort/long_cor_df.csv", row.names = F)
write.csv(short_cor_df, "/home/yoyerush/yo/methylome_pipeline/TEs_methylation_ave_cor_LongShort/short_cor_df.csv", row.names = F)



#
#
#
TE_ave_df
long_ave_df
short_ave_df
#
#
#
TE_cor_df
long_cor_df
short_cor_df