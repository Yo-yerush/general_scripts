getwd()
library(UniprotR)

uni.acce = read.csv("uniprot-pomegranate.txt", sep = "\t", header = T)

ProteinAccList = uni.acce$Cross.reference..RefSeq.
ProteinAccList = ProteinAccList[!ProteinAccList == ""]
write.table(ProteinAccList, "pome.prot.refseq.txt", row.names = F)

################
pro.acc = read.csv("pome.prot.refseq.to.uni.txt", sep = "\t", header = F)
from.pro.to.uni = ConvertID(pro.acc$V1 , ID_from = "P_REFSEQ_AC" , ID_to = "ACC", directorypath = NULL)
write.csv(from.pro.to.uni, "protein.uni.ID.csv", sep = "\t", row.names = F)

