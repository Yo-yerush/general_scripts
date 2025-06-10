setwd("P:/yonatan/RNAseq_yonatan_2021/DESeq2/uniprot/make.description.and.GO")
xm2xp = read.csv("xm.to.xp.csv", sep = "\t",header = F)
names(xm2xp) = c("transcript_id","protein_id")
uni = read.csv("uni.by.XP.csv")
uni.for.merge = uni[,c(5,4,2,3)]

uni.with.transcript = merge(xm2xp,uni.for.merge,by = "protein_id", all.x = T)

transcript.id.file = read.csv("transcript_id.and.description.csv", sep = '\t', header = F)
names(transcript.id.file) = c("transcript_id", "description")

transcript_id.and.description.with.GO.pre = merge(transcript.id.file, uni.with.transcript, all.x = T)
transcript_id.and.description.with.GO = transcript_id.and.description.with.GO.pre[-c(1:6757),c(1,3,2,4:6)]

write.csv(transcript_id.and.description.with.GO, "transcript_id.and.description.with.GO.csv", row.names = F)
