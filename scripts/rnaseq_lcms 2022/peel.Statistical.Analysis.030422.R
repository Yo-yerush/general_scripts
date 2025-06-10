library("MetaboAnalystR")

setwd("C:/Users/yonye/Downloads/yooooo/")

# PID of current job: 1098017

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "lc.T.for.metaboanalyst.csv", "colu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckData(mSet)
mSet<-FilterVariable(mSet, "none", "F", 25)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA) 
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.1, TRUE, "raw")
mSet<-PlotVolcano(mSet, "volcano_0_",1, 0, "png", 72, width=NA)
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)
mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)
mSet<-PlotTT(mSet, "tt_0_", "png", 72, width=NA)
#mSet<-PlotCorrHeatMap(mSet, "corr_0_", "png", 72, width=NA, "col", "pearson", "bwm", "overview", F, F, "0")
mSet<-PlotCorrHeatMap(mSet, "corr_1_", "png", 72, width=NA, "col", "pearson", "bwm", "detail", F, F, "0")
mSet<-PCA.Anal(mSet)
#mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
#mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_2", "png", 72, width=NA, 1,2,0.95,0,0)
#mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
#mSet<-PlotPCA3DLoading(mSet, "pca_loading3d_0_", "json", 1,2,3)
mSet<-PlotHeatMap(mSet, "heatmap_0_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "overview", T, T, NULL, T, F)
mSet<-PlotHeatMap(mSet, "heatmap_1_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", T, T, NULL, T, F)
mSet<-PreparePDFReport(mSet, "guest2401597091084768147")

file.remove(list.files()[grep("qs",list.files())])





mSet<-FeatureCorrelation(mSet, "pearson", "Pyrogallol")
mSet<-PlotCorr(mSet, "ptn_1_", "feature", "png", 72, width=NA)
mSet<-FeatureCorrelation(mSet, "pearson", "Punicalagin - STD")
mSet<-PlotCorr(mSet, "ptn_2_", "feature", "png", 72, width=NA)
mSet<-FeatureCorrelation(mSet, "pearson", "Gallic acid - STD")
mSet<-PlotCorr(mSet, "ptn_3_", "feature", "png", 72, width=NA)
mSet<-FeatureCorrelation(mSet, "pearson", "Ellagic acid - STD")
mSet<-PlotCorr(mSet, "ptn_4_", "feature", "png", 72, width=NA)
mSet<-SaveTransformedData(mSet)


