library(DMRcaller)
data(methylationDataList)
data(GEs)


#ch1 = meth_var1[meth_var1@seqnames == "NC_003071.7"]
#regions <- GRanges(seqnames = Rle("NC_003071.7"), ranges = IRanges(ch1@ranges@start[1],ch1@ranges@start[length(ch1@ranges@start)]))
regions <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(1,1E6))
# compute low resolution profile in 10 Kb windows
profileCGWT <- computeMethylationProfile(methylationDataList[["WT"]],
                                         regions,
                                         windowSize = 10000,
                                         context = "CHG")
profileCGMet13 <- computeMethylationProfile(methylationDataList[["met1-3"]],
                                            regions,
                                            windowSize = 10000,
                                            context = "CHG")
methylationProfiles <- GRangesList("WT" = profileCGWT, "met1-3" = profileCGMet13)
#plotMethylationProfile(profileCHGWT)
if (F) {
  profile_var1 <- computeMethylationProfile(meth_var1,
                                            regions,
                                            windowSize = 100000,
                                            context = "CHG")
  profile_var2 <- computeMethylationProfile(meth_var2,
                                            regions,
                                            windowSize = 100000,
                                            context = "CHG")
  methylationProfiles <- GRangesList("var1" = profile_var1, "var2" = profile_var2)
}

pch = c(1, 0, 16, 2, 15, 17)
lty = c(4,1, 3, 2, 6, 5)
col <- rainbow(length(methylationProfiles))
ymax <- 0.7

pos <- (start(methylationProfiles[[1]]) + end(methylationProfiles[[1]]))/2
plot(pos, methylationProfiles[[1]]$Proportion, type = "o", 
     ylim = c(0, ymax), xlab = "genomic coordinate", ylab = "methylation", 
     col = col[1], yaxt = "n", pch = pch[1], lty = lty[1], 
     main = "yonatan")
if (length(methylationProfiles) > 1) {
  for (i in 2:length(methylationProfiles)) {
    lines(pos, methylationProfiles[[i]]$Proportion, 
          type = "o", col = col[i], lty = lty[i], pch = pch[i])
  }
}
legend("topright", legend = names(methylationProfiles), 
       lty = lty[1:length(methylationProfiles)], col = col[1:length(methylationProfiles)], 
       pch = pch[1:length(methylationProfiles)], bty = "n")

  #mtext(labels[1], line = 0.7, adj = 0, cex = 1.4)


  axis(2, c(0, signif(ymax/2, 1), signif(ymax, 1)))



chr_local <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(5E5,6E5))

DMRsBinsCG <- computeDMRs(methylationDataList[["WT"]],
methylationDataList[["met1-3"]],
regions = chr_local,
context = "CG",
method = "bins",
binSize = 100,
test = "score",
pValueThreshold = 0.01,
minCytosinesCount = 4,
minProportionDifference = 0.4,
minGap = 200,
minSize = 50,
minReadsPerCytosine = 4,
cores = 1)

DMRsCGList <- list("bins" = DMRsBinsCG)


plotLocalMethylationProfile(methylationDataList[["WT"]],
methylationDataList[["met1-3"]],
GRanges(seqnames = Rle("Chr3"), ranges = IRanges(510000,530000)),
DMRsCGList,
conditionsNames = c("WT", "met1-3"),
GEs,
windowSize = 300,
main="CG methylation")
