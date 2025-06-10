path.to.all.gene.file = read.csv("C:/Users/yonye/OneDrive/שולחן העבודה/מעבדה/RNAseq_yonatan_2021/DESeq2/clone/statistics/BanG_VS_DanA/clone.all.BanG_VS_DanA.DE.csv")
metabo.cor.file = read.csv("C:/Users/yonye/OneDrive/שולחן העבודה/קבצים חדשים בית/250522/PatternHunter.files/results/Punicalagin B/Punicalagin B_correlation_description_table.csv")

path.to.all.gene.file = path.to.all.gene.file[,c(1,3,4)]

new.merge.file = merge.data.frame(metabo.cor.file, path.to.all.gene.file, by = "transcript_id", all.x = T)
new.merge.file = new.merge.file[,c(1,2,4,5,3)]
names(new.merge.file)[3:4] = c("log2fc","padj" )

new.merge.file = new.merge.file[order(new.merge.file$correlation, decreasing = T),]
