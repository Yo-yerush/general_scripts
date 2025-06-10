sdh_genes_cor = function(feature.name) {

ori_path = getwd()

path = "P:/yonatan/RNAseq_yonatan_2021/SDH_variants_RNAseq_yonatan/SDH_gene_puttern_hunter/" # path for output files
setwd(path)

file.name = "SDH_gene_all_high_VS_low._for_metabo_DE.csv"

#####################################################################

dir.create(paste(path, feature.name, sep = ""))
path.to.save.files = paste0(path, feature.name)


#remove.packages("tibble")
#install.packages("tibble")

library("MetaboAnalystR")
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, file.name, "colu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet)
#mSet<-SanityCheckData(mSet)
#mSet<-FilterVariable(mSet, "none", "F", 25)
mSet<-GetGroupNames(mSet, "")
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
file.remove(list.files()[grep(".qs",list.files())])

setwd(path.to.save.files)
mSet<-FeatureCorrelation(mSet, "pearson", feature.name)
mSet<-PlotCorr(mSet, paste(feature.name,"_correlation_top25_", sep = ""), "feature", "png", 150, width=NA)

corr.file = read.csv("correlation_feature.csv")
corr.file = corr.file[order(corr.file[,2], decreasing = T),]

###########
description.table = read.csv("P:/yonatan/RNAseq_yonatan_2021/lcms/PatternHunter/PatternHunter.files/transcript_id.and.description.with.GO.csv")
names(corr.file)[1] = c("transcript_id")
description.table = description.table[,c(1,3,2,4:6)]
merge.files = merge.data.frame(corr.file, description.table,by = "transcript_id", all.x = T)
merge.files = merge.files[order(-merge.files[,2]),]
###########

file.remove("correlation_feature.csv")
write.csv(merge.files, paste(feature.name,"_correlation_description_table.csv", sep = ""), row.names = F)

ori_path
rm(list=ls())
}
