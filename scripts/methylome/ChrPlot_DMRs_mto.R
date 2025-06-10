library(GenomicRanges)
res_dir = "P:/yonatan/methionine/methylome_23/BSseq_results_191223"
output_name = "mto"

for (context in c("CG","CHG","CHH")) {
  
  DMRs_mto1.0 = makeGRangesFromDataFrame(read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/mto1_vs_wt/DMRs_",context,"_mto1_vs_wt.csv")), keep.extra.columns = T)
  DMRs_mto3.0 = makeGRangesFromDataFrame(read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/mto3_vs_wt/DMRs_",context,"_mto3_vs_wt.csv")), keep.extra.columns = T)
  
  edit_dmrs <- function(df.f) {
    df.f$plog = -log10(df.f$pValue)
    df.f$plog[which(df.f$plog == Inf)] = 300
    df.f$plog[which(df.f$plog > 300)] = 300
    df.f$plog[df.f$regionType == "loss"] = -df.f$plog[df.f$regionType == "loss"]
    return(df.f)
  }

  DMRs_mto1.1 = edit_dmrs(DMRs_mto1.0)
  DMRs_mto3.1 = edit_dmrs(DMRs_mto3.0)
  
  p.max = 300
  p.min = -300
  ##################################################
  
  rename_seq <- function(gr_obj) {
    #### change seqnames function
    gr_obj_trimmed = gr_obj
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003070.9","1",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003071.7","2",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003074.8","3",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003075.7","4",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003076.8","5",seqlevels(gr_obj_trimmed)))
    return(gr_obj_trimmed)
  }
  
  DMRs_mto1 = rename_seq(DMRs_mto1.1)
  DMRs_mto3 = rename_seq(DMRs_mto3.1)
  
  chromosome_plot <- function(DMRs_trimmed.1,DMRs_trimmed.2,chr.n,ymax,ymin,pch) {
    chr.vec.1 = DMRs_trimmed.1[DMRs_trimmed.1@seqnames == chr.n]
    chr.vec.2 = DMRs_trimmed.2[DMRs_trimmed.2@seqnames == chr.n]
    
    # Subset ranges where 'plog' is between -10 and 10
    t.hold = 3#-log10(0.05)
    chr.vec.3 = c(chr.vec.1[chr.vec.1$plog >= -t.hold & chr.vec.1$plog <= t.hold],
                  chr.vec.2[chr.vec.2$plog >= -t.hold & chr.vec.2$plog <= t.hold])
    
    chr.vec.1 = c(chr.vec.1[chr.vec.1$plog <= -t.hold], chr.vec.1[chr.vec.1$plog >= t.hold])
    chr.vec.2 = c(chr.vec.2[chr.vec.2$plog <= -t.hold], chr.vec.2[chr.vec.2$plog >= t.hold])
    
    methylationProfiles <- GRangesList("var1" = chr.vec.1, "var2" = chr.vec.2, "var3" = chr.vec.3)
    
    col = c("red4","#0072B2","black")
    #pch = 20
    lty = 0
    #col_legend = col[1:length(methylationProfiles)]
    #lty_legend = lty[1:length(methylationProfiles)]
    
    pos.1 = (start(methylationProfiles[[1]]) + end(methylationProfiles[[1]]))/2
    pos.2 = (start(methylationProfiles[[2]]) + end(methylationProfiles[[2]]))/2
    pos.3 = (start(methylationProfiles[[3]]) + end(methylationProfiles[[3]]))/2
    
    plot(pos.1, methylationProfiles[[1]]$plog, type = "o",
         ylim = c(ymin, ymax), 
         col = col[1], pch = pch, lty = lty, 
         yaxt = "n", xaxt = "n", 
         #axes = F,
         main = ""
    )
    points(pos.2, methylationProfiles[[2]]$plog, col = col[2], pch = pch)
    points(pos.3, methylationProfiles[[3]]$plog, col = col[3], pch = pch)
    mtext(paste0("Chr ",chr.n), side = 1, line = 0, cex = 1)#adj = 0,
  }
  
  svg(paste0(res_dir,"/ChrPlot_DMRs_",context,"_",output_name,".svg"), width = 8, height = 2, family = "serif")
  
  #left axis
  par(mar = c(1,4,2,0))
  par(fig=c(0,2,0,10)/10)
  plot(runif(10),runif(10),xlim=c(0, 0.01), 
       ylim=c(p.min,p.max),axes=FALSE,type="n",ylab=paste(context," DMRs (-log[p])"),xlab="")
  axis(2, c(p.min, 0, p.max), labels = c(-p.min, 0, p.max))
  par(new=T)
  
  # main ChrPlot
  i=1
  i.0=1.15
  u = (10-i.0)/6
  for (chr_number in 1:5) {
    par(mar = c(1,0,2,0))
    par(fig=c(i,i+u,0,10)/10)
    chromosome_plot(DMRs_mto1,DMRs_mto3,chr_number,p.max,p.min,"•")
    par(new=T)
    i=i+u
  }
  
  # right axis
  par(mar = c(1,0,2,2))
  par(fig=c(u*4+1-0.025,u*4+3-0.025,0,10)/10)
  axis(4, c(p.min, 0, p.max), labels = F, )
  par(new=T)
  
  # right axis - gain or loss
  par(mar = c(1,0,2,2))
  par(fig=c(6.6,8.6,0,10)/10) # c(u*4+1,u*4+3,0,10)/10
  axis(4, c(p.min/2, p.max/2), labels = c("Loss", "Gain"), tick = F)
  par(new=T)
  
  # legend
  par(mar = c(1,2,2,0))
  par(fig=c(8,10,0,10)/10) # c(u*4+3,10,0,10)/10)
  legend("topright", legend = c("mot1 vs wt","mot3 vs wt"), 
         lty = 0, col = c("red4","#0072B2"), 
         pch = 16, bty = "n")
  dev.off()
  ##################################################
}

