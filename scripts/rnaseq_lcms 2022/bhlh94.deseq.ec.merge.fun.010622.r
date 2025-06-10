ec_deseq_merge = function(peel.file){

setwd("P:/yonatan/RNAseq_yonatan_2021/for.enzym/ec_2.3.1/")

pro.list = read.csv("pom_uniprot_ec_2.3.1_300522.csv", header = F)
names(pro.list) = "protein_id"
all.gene.exp.file = read.csv("P:/yonatan/RNAseq_yonatan_2021/DESeq2/transgenic/statistics/bHLH94-like_VS_EV/transgenic.all.bHLH94-like_VS_EV.DE.csv")[,-c(6,8)]

only.EC.file = merge.data.frame(all.gene.exp.file, pro.list, by = "protein_id", all.y = T)
only.EC.file = only.EC.file[order(only.EC.file$padj),]
#only.EC.file = only.EC.file[only.EC.file$padj <= 0.05,]

write.csv(only.EC.file, paste("peel.deseq.ec.merge.files/bHLH94-like_VS_EV_ec.2.3.1_merge.csv", sep = ""), row.names = F)
}