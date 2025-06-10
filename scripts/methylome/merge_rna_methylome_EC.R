library(dplyr)
library(writexl)
library(stringr)
library(RColorBrewer)
library(rtracklayer)
library(circlize)


des_file = read.csv("P:/TAIR10.1/AT_description_file_260123.csv")

# find ec 2.1.1.-
all_interesting_tairs = unique(des_file[grep("^2\\.1\\.1\\.", des_file$EC.number) ,"locus_tag"])


xl_list = list()

for (treatment in c("mto1","mto3")) {
  #ann.m = ifelse(ann == "genes", "Genes","Promoters")
  ############# methylome and RNAseq file filtered by corr files #############
  #corr_file = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/corr_with_methylations/",treatment,"/",context,"/",ann,".corr.",context,".",treatment,".csv"))
  #corr_file = corr_file[corr_file$pval < 0.05, c("transcript_id","cor")]
  
  
  RNA_file = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/met23/",treatment,"_vs_wt/all.transcripts.",treatment,"_vs_wt.DE.csv"))
  #RNA_file = RNA_file[RNA_file$padj < 0.05, c("transcript_id","log2FoldChange")]# %>% na.omit()
  RNA_file = RNA_file[, c("transcript_id","locus_tag","log2FoldChange","padj")]
  names(RNA_file)[3:4] = c("RNA_log2FC","RNA_padj")
  RNA_file$RNA_log2FC = round(RNA_file$RNA_log2FC, 3)
  RNA_file$RNA_padj = as.numeric(sprintf("%.3e", RNA_file$RNA_padj))
  
  #merge1 = merge.data.frame(corr_file,RNA_file, by = "transcript_id", all.y = T)
  meth_fun = function(treatment.f, context.f, ann.f) {
    meth_file.f = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment.f,"_vs_wt/genome_annotation/",context.f,"/",ann.f,"_",context.f,"_genom_annotations.csv"))
    meth_file.f$log2FoldChange = round(log2(meth_file.f$proportion2 / meth_file.f$proportion1), 3)
    #meth_file.f = meth_file.f[,c("locus_tag","seqnames","start","end","log2FoldChange")]
    meth_file.f = meth_file.f[,c("locus_tag","log2FoldChange")]
    names(meth_file.f)[2] = paste0(context.f,"_",ann.f)
    return(meth_file.f)
  }
  
  
  meth_list = list()
  meth_list[["CG_genes"]] = meth_fun(treatment,"CG","Genes")
  meth_list[["CHG_genes"]] = meth_fun(treatment,"CHG","Genes")
  meth_list[["CHH_genes"]] = meth_fun(treatment,"CHH","Genes")
  meth_list[["CG_promoters"]] = meth_fun(treatment,"CG","Promoters")
  meth_list[["CHG_promoters"]] = meth_fun(treatment,"CHG","Promoters")
  meth_list[["CHH_promoters"]] = meth_fun(treatment,"CHH","Promoters")
  
  merged_meth = meth_list[[1]]
  for (i in 2:length(meth_list)) {
    merged_meth = merge.data.frame(merged_meth, meth_list[[i]], by = "locus_tag", all = TRUE)
  }
  
  merge2 = merge.data.frame(RNA_file, merged_meth, by = "locus_tag", all.x = T)

  filtered_merge2 = merge2[grepl(paste(all_interesting_tairs, collapse = "|"), merge2$locus_tag),]

  merge3 = merge.data.frame(filtered_merge2, des_file[,-1], by = "transcript_id", all.x = T)
  merge3 <- merge3[order(merge3$RNA_padj), ]
  #merge3 <- merge3[order(merge3$cor, decreasing = T), ]
  
  xl_list[[treatment]] = merge3
}

write_xlsx(xl_list, "P:/yonatan/methionine/merged_results_mtos_EC_250123.xlsx")
#