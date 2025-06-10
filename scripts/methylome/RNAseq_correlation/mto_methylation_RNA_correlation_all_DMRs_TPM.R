library(plyr)
library(dplyr)
library(corrplot)
library(ggpubr)
library(ggplot2)
library(lmerTest)
library(tidyr)

for (treatment in c("mto1","mto3")) {
  for (context in c("CG","CHG","CHH")) {
    for (annotation_type in c("genes","promoters")) {
      
      path_2_save.0 = "P:/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_TPM/"
      path_2_save.1 = paste0(path_2_save.0,treatment,"/")
      path_2_save.2 = paste0(path_2_save.1,context,"/")
      dir.create(path_2_save.1)
      dir.create(path_2_save.2)
      
      if(treatment == "mto1") {
        mto_rnaseq_names = c("met14","met15","met16")
      } else {
        mto_rnaseq_names = c("met17","met18","met19")
      }
      wt_rnaseq_names = c("met20","met22")
      
      #Load the methylation files of the promoters of interestL
      meth_matrix = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/average.meth.genes.levels/",treatment,"_vs_wt/",context,"/meth.",annotation_type,".",context,".",treatment,"_vs_wt.csv"))
      
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
      
      #RNA.0 = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/met23/",treatment,"_vs_wt/norm.",treatment,"_vs_wt.DE.csv"))

      #### Gene expression ####
      refseq_2_tair = read.csv("P:/TAIR10.1/TAIR_2_refseq_transcript.csv")[,1:2]
      RNA = merge.data.frame(refseq_2_tair, RNA.0, by = "transcript_id")
      # edit table and names
      mto_names_pos = grep(paste(mto_rnaseq_names, collapse = "|"), names(RNA))
      names(RNA)[mto_names_pos] = paste0(treatment,".",1:3,"_RNA")
      
      wt_names_pos = grep(paste(wt_rnaseq_names, collapse = "|"), names(RNA))
      names(RNA)[wt_names_pos] = paste0("wt.",1:2,"_RNA")
      
      RNA = RNA[,c(1,2,wt_names_pos,mto_names_pos)]
      
      
      #### Methylation ####
      meth_matrix = meth_matrix[,-c(1:5)]
      names(meth_matrix)[-(1:2)] = paste0(names(meth_matrix)[-(1:2)], "_meth")
      
      #Merge methylation and RNA dataframes
      meth_matrix_RNA = merge(meth_matrix, RNA, by="locus_tag", all=FALSE)
      meth_matrix_RNA <- meth_matrix_RNA %>%
        select(locus_tag, transcript_id, everything())
      
      meth_new = meth_matrix_RNA[,c(grep("locus_tag|gene",names(meth_matrix_RNA)), grep("meth",names(meth_matrix_RNA)))]
      RNA_new = meth_matrix_RNA[,c(grep("locus_tag|transcript_id",names(meth_matrix_RNA)), grep("RNA",names(meth_matrix_RNA)))]
      
      if (all(meth_new$locus_tag == RNA_new$locus_tag)) {
        merge_new = cbind(RNA_new[,-1], meth_new[,-(1:2)]) # for cor_df to save
      }
      
      names(meth_new) = gsub("_meth","",names(meth_new))
      names(RNA_new) = gsub("_RNA","",names(RNA_new))
      
      #move the columns/rows:
      fin_meth = meth_new %>%
        pivot_longer(!c(locus_tag,gene), names_to = "Sample", values_to = "Meth")
      
      fin_RNA = RNA_new %>%
        pivot_longer(!c(locus_tag,transcript_id), names_to = "Sample", values_to = "TPM")
      
      # check if both data frames are similar (for cbind)
      if (
        all(fin_meth$Sample == fin_RNA$Sample) & all(fin_meth$locus_tag == fin_RNA$locus_tag)
      ) {
        
        fin_meth_RNA = cbind(fin_meth, fin_RNA[,c("transcript_id","TPM")])
        fin_meth_RNA <- fin_meth_RNA %>%
          select(locus_tag, transcript_id, everything())
        fin_meth_RNA$Meth[is.na(fin_meth_RNA$Meth)] = 0 # i dont know if its the right way
        
        
        ## Run the correlation analysis and correct for multiple comparisons (FDR < 0.05)
        cor_df = fin_meth_RNA %>%
          group_by(transcript_id) %>% # check which column to group: 'locus_tag' or 'transcript_id'
          summarise(cor = stats::cor.test(Meth, TPM, method="spearman", exact=F)$estimate,
                    pval = stats::cor.test(Meth, TPM, method="spearman", exact=F)$p.value) %>%
          ungroup() %>% na.omit()
        cor_df = merge.data.frame(fin_meth_RNA[,1:3], cor_df, by = "transcript_id")
        cor_df <- cor_df %>% 
          distinct(transcript_id, .keep_all = TRUE)
        
        #cor_df = cor_df[(cor_df[,"pval"] < 0.05),]
        #cor_df = na.omit(cor_df)
        cor_df$padj = p.adjust(cor_df$pval, method = "fdr", n = nrow(cor_df))
        
        cor_df_2_save = merge.data.frame(cor_df, merge_new, by = "transcript_id") %>%
          arrange(pval)
        
        message(nrow(cor_df_2_save[(cor_df_2_save[,"pval"] < 0.05),]), " significant correlated ",annotation_type," (",context," - ",treatment,")")
        write.csv(cor_df_2_save, paste0(path_2_save.2,annotation_type,".corr.",context,".",treatment,".csv"), row.names = F)
      }
      
    }
  }
}