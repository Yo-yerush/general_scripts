setwd("P:/yonatan/RNAseq_yonatan_2021/DESeq2/uniprot")

uni.acce = read.csv("uniprot-pomegranate.txt", sep = "\t", header = T)
uni.pro = uni.acce[,c(7,6,4,5)]
names(uni.pro) = c("protein_id", "biological.process", "molecular.function", "cellular.component")
uni.pro = uni.pro[!uni.pro$protein_id == "",]


#uni.pro2 = uni.pro[grep(".1;X" ,uni.pro$protein_id),]
#View(uni.pro2)



trans.prot.1 = read.csv("transcript.protein.description.csv", sep = "\t", header = F)
trans.prot = trans.prot.1[,c(3,1,2)]
names(trans.prot) = c("protein_id","transcript_id", "description")

id.file = merge(trans.prot, uni.pro, by = "protein_id", all = T)
View(id.file)


#### need to edit protein_id from uniprot. ( ; )
