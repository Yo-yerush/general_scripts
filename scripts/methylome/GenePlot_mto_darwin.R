methylationData_mto1.1 <- readBismark("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S9/methylation_extractor/S9_R1_bismark_bt2_pe.CX_report.txt")
methylationData_mto1.2 <- readBismark("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S10/methylation_extractor/S10_R1_bismark_bt2_pe.CX_report.txt")
methylationData_mto1.3 <- readBismark("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S11/methylation_extractor/S11_R1_bismark_bt2_pe.CX_report.txt")

methylationData_mto3.1 <- readBismark("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S12/methylation_extractor/S12_R1_bismark_bt2_pe.CX_report.txt")
methylationData_mto3.2 <- readBismark("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S13/methylation_extractor/S13_R1_bismark_bt2_pe.CX_report.txt")
methylationData_mto3.3 <- readBismark("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S14/methylation_extractor/S14_R1_bismark_bt2_pe.CX_report.txt")

methylationData_wt.1 <- readBismark("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S18/methylation_extractor/S18_R1_bismark_bt2_pe.CX_report.txt")
methylationData_wt.2 <- readBismark("/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S19/methylation_extractor/S19_R1_bismark_bt2_pe.CX_report.txt")

methylationPool_mto1 <- poolMethylationDatasets(GRangesList(mto1.1 = methylationData_mto1.1,
                                                            mto1.2 = methylationData_mto1.2,
                                                            mto1.3 = methylationData_mto1.3))

methylationPool_mto3 <- poolMethylationDatasets(GRangesList(mto3.1 = methylationData_mto3.1,
                                                            mto3.2 = methylationData_mto3.2,
                                                            mto3.3 = methylationData_mto3.3))

methylationPool_wt <- poolMethylationDatasets(GRangesList(wt.1 = methylationData_wt.1,
                                                          wt.2 = methylationData_wt.2))
##################################################
#### change seqnames function
rename_seq <- function(gr_obj) {
  gr_obj_trimmed = gr_obj
  gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003070.9","1",seqlevels(gr_obj_trimmed)))
  gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003071.7","2",seqlevels(gr_obj_trimmed)))
  gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003074.8","3",seqlevels(gr_obj_trimmed)))
  gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003075.7","4",seqlevels(gr_obj_trimmed)))
  gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003076.8","5",seqlevels(gr_obj_trimmed)))
  gr_obj_trimmed = gr_obj_trimmed[!gr_obj_trimmed@seqnames == "NC_037304.1"]
  gr_obj_trimmed = gr_obj_trimmed[!gr_obj_trimmed@seqnames == "NC_000932.1"]
  return(gr_obj_trimmed)
}
##################################################
chromosome_plot <- function(meth_var1_trimmed,meth_var2_trimmed,meth_var3_trimmed,#var1_name,var2_name,
                            chr.n,cntx,ymax,ymin) {
  
  chr.vec = meth_var1_trimmed[meth_var1_trimmed@seqnames == chr.n]
  regions <- GRanges(seqnames = Rle(chr.n), ranges = IRanges(chr.vec@ranges@start[1],chr.vec@ranges@start[length(chr.vec@ranges@start)]))
  
  profile_var1 <- computeMethylationProfile(meth_var1_trimmed,
                                            regions,
                                            windowSize = 150000,
                                            context = cntx)
  
  profile_var2 <- computeMethylationProfile(meth_var2_trimmed,
                                            regions,
                                            windowSize = 150000,
                                            context = cntx)
  profile_var3 <- computeMethylationProfile(meth_var3_trimmed,
                                          regions,
                                          windowSize = 150000,
                                          context = cntx)
  
  methylationProfiles_var1 = profile_var1[,0]
  methylationProfiles_var1$Proportion = profile_var1$Proportion - profile_var3$Proportion
  
  methylationProfiles_var2 = profile_var2[,0]
  methylationProfiles_var2$Proportion = profile_var2$Proportion - profile_var3$Proportion
  
  methylationProfiles <- GRangesList("var1" = methylationProfiles_var1, "var2" = methylationProfiles_var2)
  
  #legend_text = c(var1_name,var2_name,var3_name)
  col = c("red4","#0072B2","#e37320")
  pch = rep(26,length(methylationProfiles))
  lty = rep(1,length(methylationProfiles))
  col_legend = col[1:length(methylationProfiles)]
  lty_legend = lty[1:length(methylationProfiles)]
  
  pos = (start(methylationProfiles[[1]]) + end(methylationProfiles[[1]]))/2
  plot(pos, methylationProfiles[[1]]$Proportion, type = "o", 
       ylim = c(ymin, ymax), 
       col = col[1], pch = pch[1], lty = lty[1], 
       yaxt = "n", xaxt = "n", 
       #axes = F,
       #main = paste0("Methylation in ",cntx," context")
       
       main = ""#paste0("Chromosome ",chr.n)
  )
  if (length(methylationProfiles) > 1) {
    for (i in 2:length(methylationProfiles)) {
      lines(pos, methylationProfiles[[i]]$Proportion, 
            type = "o", col = col[i], lty = lty[i], pch = pch[i])
    }
  }
  lines(pos, rep(0, length(methylationProfiles[[1]])), 
        type = "o", col = "#141414", lty = 1, pch = 26)
  #legend("topright", legend = legend_text, 
  #       lty = lty_legend, col = col_legend, 
  #       pch = pch[1:length(methylationProfiles)], bty = "n")
  mtext(paste0("Chr ",chr.n), side = 1, line = 0, cex = 1)#adj = 0,
  
}
##################################################
meth_mto1_trimmed = rename_seq(methylationPool_mto1)
meth_mto3_trimmed = rename_seq(methylationPool_mto3)
meth_wt_trimmed = rename_seq(methylationPool_wt)
##################################################
for (cntx in c("CG","CHG","CHH")) {

  if (cntx == "CG") {
    ymax = 0.1
    ymin = -0.1
  } else if (cntx == "CHG") {
    ymax = 0.1
    ymin = -0.1
  } else if (cntx == "CHH") {
    ymax = 0.03
    ymin = -0.03
  }
  
  #### difference plot
  svg(paste0("ChrPlot_difference_",cntx,"_mto",".svg"), width = 8, height = 2, family = "serif")
  par(mar = c(1,4,2,0))
  par(fig=c(0,2,0,10)/10)
  plot(runif(10),runif(10),xlim=c(0, 0.01), 
       ylim=c(ymin,ymax),axes=FALSE,type="n",ylab=paste(cntx," methylation (Δ)"),xlab="")
  axis(2, c(ymin, 0, ymax), lty = 1)
  par(new=T)
  
  i=1
  u = (10-i)/5
  for (chr_number in 1:5) {
    par(mar = c(1,0,2,0))
    par(fig=c(i,i+u,0,10)/10)
    chromosome_plot(meth_mto1_trimmed,meth_mto3_trimmed,meth_wt_trimmed,#var1_name,var2_name,
                    chr_number,cntx,ymax,ymin)
    par(new=T)
    i=i+u
  }
  dev.off()
}