# 13/04/2022
# Author: Iris Sammarco
# Description: prepare the input files to calculate correlation between promoter methylation and gene expression

library(plyr)
library(dplyr)
library(corrplot)
library(ggpubr)
library(ggplot2)
library(lmerTest)
library(tidyr)

#set your working directory here
setwd("")

#### Gene expression ####
RNA = read.table("FPKM_all_samples_noNAs_IDs.bedgraph", sep= "\t", dec = ".", header = TRUE) # RNA input file, it looks like this:
#gene_ID	sample1	sample2
#FvH4_1g00020	0.026454434	0.026777595
#FvH4_1g00021	0.04		0.1
#the numbers are FPKM values per each gene of interest per sample
RNA$gene_ID = as.character(RNA$gene_ID)

#### CpG ####
## Methylation
#Load the methylation files of the promoters of interestL
CpG_meth_matrix = read.table("CpG_G_7vs7_avg_meth_env_all_DMR_promoters_IDs_direct.bedgraph", sep= "\t", dec = ".", header = TRUE) #this file looks like this:
#chrom	start	end	sample1	sample2	gene_ID
#Fvb1	16158	16280	0.584699453551913	0.597315436241611	FvH4_1g00020
#Fvb1	17000	17200	0.7	0.6	FvH4_1g00021
#the numbers are methylation frequencies each gene promoter of interest per sample

#Merge methylation and RNA dataframes
CpG_meth_matrix_RNA = merge(CpG_meth_matrix, RNA, by="gene_ID", all=FALSE)
#rename the first cols:
colnames(CpG_meth_matrix_RNA)[colnames(CpG_meth_matrix_RNA) == "chr.x"] = "chr"
colnames(CpG_meth_matrix_RNA)[colnames(CpG_meth_matrix_RNA) == "start.x"] = "start"
colnames(CpG_meth_matrix_RNA)[colnames(CpG_meth_matrix_RNA) == "end.x"] = "end"
write.csv(CpG_meth_matrix_RNA, "CpG_G_7vs7_avg_meth_direct_env_DMR_promoters_FPKM.csv", row.names = FALSE) #you can check how the file looks like

## I want to create a df that will allow me to test for correlation between methylation and expression for each gene and for all the samples, so I want a df that looks like this:
#gene_ID Meth FPKM Sample
#FvH4_1g00020 0.6 0.2 sample1
#FvH4_1g00020 0.8 0.1 sample2

# If you have different samples for methylation and RNA data, select only the methylation data for the samples for which you have both meth and RNA data, in this way I know that for the same number of rows, I'm selecting the same genes for both RNA and meth:
## Methylation df:
CpG_meth_matrix_sub_samples= CpG_meth_matrix_RNA[1:67]
#remove chr start end cols:
CpG_meth_matrix_sub_samples= CpG_meth_matrix_RNA[, -c(2:4)]
#move the columns/rows:
fin_meth = CpG_meth_matrix_sub_samples %>%
  pivot_longer(!gene_ID, names_to = "Sample", values_to = "Meth")

## RNA df:
CpG_RNA_matrix_sub_samples=CpG_meth_matrix_RNA[, -c(2:70)]
#move the columns/rows:
fin_RNA = CpG_RNA_matrix_sub_samples %>%
  pivot_longer(!gene_ID, names_to = "Sample", values_to = "FPKM")

# Merge the meth and rna dfs (the genes are in the same order):
fin_meth_RNA = cbind(fin_meth, fin_RNA)
#check whether the gene_ID from RNA and meth are identical:
all(fin_meth_RNA$V1 == fin_meth_RNA$V4) #TRUE: they are identical!
#drop duplicated columns:
fin_meth_RNA =fin_meth_RNA[,-c(4,5)]
write.csv(fin_meth_RNA, "CpG_G_7vs7_avg_meth_direct_env_DMR_promoters_FPKM_df.csv", row.names = FALSE)

## Run the correlation analysis and correct for multiple comparisons (FDR < 0.05)
cor_df = fin_meth_RNA %>%
  group_by(gene_ID) %>%
  summarise(cor = stats::cor.test(Meth, FPKM, method=c("spearman"))$estimate,
            pval = stats::cor.test(Meth, FPKM, method=c("spearman"))$p.value
  ) %>%
  ungroup()
cor_df = na.omit(cor_df) #remove NAs
#select the genes significantly correlated p < 0.05 (without correcting for multiple comparisons):
cor_sign_nofdr = cor_df[(cor_df[,3] < 0.05),]
nrow(cor_sign_nofdr) #number of significant correlations
write.csv(cor_sign_nofdr, "no_FDR/CpG_G_corr_env_DMR_promoters_FPKM.csv", row.names = FALSE)
#adjust the p-values for multiple comparisons (fdr):
padj = p.adjust(cor_df$pval, method = c("fdr"), n = length(cor_df$pval))
cor_df_p = cbind(cor_df, padj)
#select only genes with FDR < 0.05
cor_sign = cor_df_p[(cor_df_p[,4] < 0.05),]
nrow(cor_sign) #number of significant correlations
write.csv(cor_sign, "CpG_G_corr_env_DMR_promoters_FPKM_FDR0.05.csv", row.names = FALSE)

# Repeat the same analysis for CHG and CHH methylation



