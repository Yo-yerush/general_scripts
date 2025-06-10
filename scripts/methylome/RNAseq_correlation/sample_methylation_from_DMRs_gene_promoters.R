library(DMRcaller)
library(parallel)
scripts_dir = "/home/yoyerush/yo/methylome_pipeline/Methylome.At/Methylome.At_scripts/"

wt_path = read.table("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_mto1.txt", sep = "\t")[1:2,2]
mto1_path = read.table("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_mto1.txt", sep = "\t")[3:5,2]
mto3_path = read.table("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_mto3.txt", sep = "\t")[3:5,2]

source(paste0(scripts_dir,"load_replicates.R"))
load_vars = mclapply(list(wt_path,mto1_path,mto3_path), function(x) load_replicates(x, 10), mc.cores = 3)

trimm_Chr <- function(gr_obj) {
  remove_seqnames = c("NC_000932.1","NC_037304.1")
  gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
  seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
  return(sort(gr_obj))
}
meth_wt = trimm_Chr(load_vars[[1]]$methylationDataReplicates)
meth_mto1 = trimm_Chr(load_vars[[2]]$methylationDataReplicates)
meth_mto3 = trimm_Chr(load_vars[[3]]$methylationDataReplicates)

meth_wt$wt.1 = meth_wt$readsM1/meth_wt$readsN1
meth_wt$wt.2 = meth_wt$readsM2/meth_wt$readsN2
meth_mto1$mto1.1 = meth_mto1$readsM1/meth_mto1$readsN1
meth_mto1$mto1.2 = meth_mto1$readsM2/meth_mto1$readsN2
meth_mto1$mto1.3 = meth_mto1$readsM3/meth_mto1$readsN3
meth_mto3$mto3.1 = meth_mto3$readsM1/meth_mto3$readsN1
meth_mto3$mto3.2 = meth_mto3$readsM2/meth_mto3$readsN2
meth_mto3$mto3.3 = meth_mto3$readsM3/meth_mto3$readsN3

meth_GRlist = list(wt = meth_wt, mto1 = meth_mto1, mto3 = meth_mto3)
############################### after it make a loop
ave.meth.genes.levels <- function(meth_ctrl, meth_trmt, treatment, cntx) {
  
  path.0 = "/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/"
  path.1 = paste0(path.0,"average.DMR.genes.levels")
  path.2 = paste0(path.1,"/",treatment)
  path.3 = paste0(path.2,"/",cntx)
  dir.create(path.1, showWarnings = F)
  dir.create(path.2, showWarnings = F)
  dir.create(path.3, showWarnings = F)
  
  meth_ctrl_cntx = meth_ctrl[which(meth_ctrl$context == cntx)]
  meth_trmt_cntx = meth_trmt[which(meth_trmt$context == cntx)]
  
  meth_ctrl_cntx = meth_ctrl_cntx[,grep("wt", names(meth_ctrl_cntx@elementMetadata))]
  meth_trmt_cntx = meth_trmt_cntx[,grep(strsplit(treatment, "_")[[1]], names(meth_trmt_cntx@elementMetadata))]
  
  keep_cols = c(2,3,4,5,6,1,8)
  DMRs_genes = makeGRangesFromDataFrame(read.csv(paste0("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/",treatment,"/genome_annotation/",cntx,"/Genes_",cntx,"_genom_annotations.csv"))[,keep_cols], keep.extra.columns = T)
  DMRs_promoters = makeGRangesFromDataFrame(read.csv(paste0("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/",treatment,"/genome_annotation/",cntx,"/Promoters_",cntx,"_genom_annotations.csv"))[,keep_cols], keep.extra.columns = T)
  
  # find overlaps
  dmrs_list = list(DMRs_genes = DMRs_genes, 
                   DMRs_promoters = DMRs_promoters)
  
  meth_list = list(meth_ctrl_cntx = meth_ctrl_cntx, 
                   meth_trmt_cntx = meth_trmt_cntx)
  
  fun4parallel <- function(ML, dmrs_list, DL) {
    m <- findOverlaps(dmrs_list[[DL]],meth_list[[ML]])
    DMRs_annotation <- dmrs_list[[DL]][queryHits(m)]
    mcols(DMRs_annotation) <- cbind.data.frame(mcols(DMRs_annotation),
                                               mcols(meth_list[[ML]][subjectHits(m)]))
    
    split_gr <- split(DMRs_annotation, DMRs_annotation$locus_tag)
    meth.0 = GRanges()
    for (i.l in 1:length(split_gr)) {
      meth.0[i.l] = split_gr[[i.l]][1]
      
      for (col.numbers in grep("\\.[1-9]",names(meth.0@elementMetadata))) {
        mcols(meth.0)[i.l, col.numbers] =  mean(mcols(split_gr[[i.l]])[,col.numbers], na.rm=TRUE)
      }
      #cat(i.l,"/",length(split_gr),"\n")
    }
    return(meth.0)
  }
  
  meth.list <- list()
  for (DL in 1:length(dmrs_list)) {
    meth.list[[DL]] <- mclapply(1:length(meth_list), function(ML) fun4parallel(ML, dmrs_list, DL), mc.cores=2)
  }
  
  meth.genes.ctrl = meth.list[[1]][[1]]
  meth.genes.trmt = meth.list[[1]][[2]]
  meth.promoters.ctrl = meth.list[[2]][[1]]
  meth.promoters.trmt = meth.list[[2]][[2]]
  
  mcols(meth.genes.ctrl) = cbind(mcols(meth.genes.ctrl),mcols(meth.genes.trmt)[,-c(1:2)])
  mcols(meth.promoters.ctrl) = cbind(mcols(meth.promoters.ctrl),mcols(meth.promoters.trmt)[,-c(1:2)])
  
  write.csv(meth.genes.ctrl, paste0(path.3,"/meth.genes.",cntx,".",treatment,".csv"), row.names = F)
  write.csv(meth.promoters.ctrl, paste0(path.3,"/meth.promoters.",cntx,".",treatment,".csv"), row.names = F)
  
  cat("\n**\ndone :",treatment," - ",cntx,"\n**\n")
}

outer_results <- mclapply(c("mto1", "mto3"), function(treatment.loop) {
  inner_results <- mclapply(c("CG", "CHG", "CHH"), function(cntx.loop) {
    ave.meth.genes.levels(meth_GRlist[["wt"]], meth_GRlist[[treatment.loop]], paste0(treatment.loop, "_vs_wt"), cntx.loop)
  }, mc.cores = 3)
  return(inner_results)
}, mc.cores = 2)
#