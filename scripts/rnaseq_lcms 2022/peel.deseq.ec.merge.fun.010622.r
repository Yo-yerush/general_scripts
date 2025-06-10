ec_deseq_merge = function(peel.file, save.files = T){

setwd("P:/yonatan/RNAseq_yonatan_2021/for.enzym/ec_2.3.1/")

pro.list = read.csv("pom_uniprot_ec_2.3.1_300522.csv", header = F)
names(pro.list) = "protein_id"
all.gene.exp.file = read.csv(paste("../../DESeq2/peel/statistics/", peel.file, "/peel.all.", peel.file, ".DE.csv", sep = ""))[,-c(6,8)]

only.EC.file = merge.data.frame(all.gene.exp.file, pro.list, by = "protein_id", all.y = T)
only.EC.file = only.EC.file[order(only.EC.file$padj),]
#only.EC.file = only.EC.file[only.EC.file$padj <= 0.05,]

if (save.files == T) {
   write.csv(only.EC.file, paste("peel.deseq.ec.merge.files/", peel.file, "_ec.2.3.1_merge.csv", sep = ""), row.names = F)
}
}