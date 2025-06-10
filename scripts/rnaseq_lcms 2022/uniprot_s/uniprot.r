library(UniprotR)

path = setwd("P:/yonatan/RNAseq_yonatan_2021/DESeq2/uniprot")

uni.acce = read.csv("uniprot-pomegranate.txt", sep = "\t", header = T)
ProteinAccList = uni.acce$Entry

from.uni.to.trans = ConvertID(ProteinAccList , ID_from = "ACC+ID" , ID_to = "REFSEQ_NT_ID", directorypath = NULL)

uni.file = from.uni.to.trans[grep("M",from.uni.to.trans$`To REFSEQ_NT_ID`),]

write.table(uni.file, "uniprot_converted.ID.pome.txt", sep = "\t", row.names = F)








