library(data.table)

setwd("P:/yonatan/RNAseq_yonatan_2021/DESeq2/uniprot/uni.files.forR")
uni.files = list.files(path = getwd(), all.files = T)
bind.file = rbindlist(lapply(uni.files[-c(1:2)], fread))
names(bind.file) = c("Entry","molecular function","cellular component","biological process","protein_id")
setwd("P:/yonatan/RNAseq_yonatan_2021/DESeq2/uniprot")
write.csv(bind.file, "uni.by.XP.csv",row.names = F)
