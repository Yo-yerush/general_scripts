library(DMRcaller)
library(parallel)

wt_path = read.table("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_mto1.txt", sep = "\t")[1:2,2]
mto1_path = read.table("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_mto1.txt", sep = "\t")[3:5,2]
mto3_path = read.table("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_mto3.txt", sep = "\t")[3:5,2]

source("/home/yoyerush/yo/methylome_pipeline/Methylome.At/scripts/load_replicates.R")
source("/home/yoyerush/yo/methylome_pipeline/Methylome.At/scripts/trimm_and_rename_seq.R")
load_vars = mclapply(list(wt_path,mto1_path,mto3_path), function(x) load_replicates(x, 10), mc.cores = 3)

meth_wt = trimm_and_rename(load_vars[[1]]$methylationDataReplicates)
meth_mto1 = trimm_and_rename(load_vars[[2]]$methylationDataReplicates)
meth_mto3 = trimm_and_rename(load_vars[[3]]$methylationDataReplicates)

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
  
  path.0 = "/home/yoyerush/yo/methylome_pipeline/"
  path.1 = paste0(path.0,"average.meth.genes.levels")
  path.2 = paste0(path.1,"/",treatment)
  path.3 = paste0(path.2,"/",cntx)
  dir.create(path.1, showWarnings = F)
  dir.create(path.2, showWarnings = F)
  dir.create(path.3, showWarnings = F)
  
  meth_ctrl_cntx = meth_ctrl[which(meth_ctrl$context == cntx)]
  meth_trmt_cntx = meth_trmt[which(meth_trmt$context == cntx)]
  
  meth_ctrl_cntx = meth_ctrl_cntx[,grep("wt", names(meth_ctrl_cntx@elementMetadata))]
  meth_trmt_cntx = meth_trmt_cntx[,grep(strsplit(treatment, "_")[[1]], names(meth_trmt_cntx@elementMetadata))]
  
  keep_cols = c(2,3,4,5,6,1,13)
  # Specify the annotation types
  annotation_types <- c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs", "TEG")

  # Read in each annotation type in a loop
  dmrs_list <- lapply(annotation_types, function(region) {
      csv_file <- paste0(
          "/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/",
          treatment, "/genome_annotation/", cntx, "/", region, "_", cntx, "_genom_annotations.csv"
      )
      makeGRangesFromDataFrame(
          read.csv(csv_file)[, keep_cols],
          keep.extra.columns = TRUE
      )
  })
  names(dmrs_list) <- annotation_types

  # This is your list of control/treatment methylation data
  meth_list <- list(
      meth_ctrl_cntx = meth_ctrl_cntx,
      meth_trmt_cntx = meth_trmt_cntx
  )

  # Define the overlap function (unchanged)
  fun4parallel <- function(ML, dmrs_list, DL) {
      m <- findOverlaps(dmrs_list[[DL]], meth_list[[ML]])
      DMRs_annotation <- dmrs_list[[DL]][queryHits(m)]
      mcols(DMRs_annotation) <- cbind.data.frame(
          mcols(DMRs_annotation),
          mcols(meth_list[[ML]][subjectHits(m)])
      )

      split_gr <- split(DMRs_annotation, DMRs_annotation$gene_id)
      meth.0 <- GRanges()

      for (i.l in seq_along(split_gr)) {
          meth.0[i.l] <- split_gr[[i.l]][1]
          for (col.numbers in grep("\\.[1-9]", names(meth.0@elementMetadata))) {
              mcols(meth.0)[i.l, col.numbers] <- mean(
                  mcols(split_gr[[i.l]])[, col.numbers],
                  na.rm = TRUE
              )
          }
          cat(i.l, "/", max(seq_along(split_gr)), "\n")
      }
      return(meth.0)
  }

  # Process all annotation types in parallel
  meth.list <- lapply(seq_along(dmrs_list), function(DL) {
      mclapply(seq_along(meth_list),
          function(ML) fun4parallel(ML, dmrs_list, DL),
          mc.cores = 2
      )
  })

  # Now 'meth.list' is a list of length = #annotation_types,
  # where each element is a length-2 list (control & treatment).
  # We can write everything in one loop:

  for (i in seq_along(annotation_types)) {
      region <- annotation_types[i]

      # Control & treatment are each a GRanges object
      region_ctrl <- meth.list[[i]][[1]]
      region_trmt <- meth.list[[i]][[2]]

      # Combine control & treatment columns (minus the first two columns, as in your original code)
      mcols(region_ctrl) <- cbind(
          mcols(region_ctrl),
          mcols(region_trmt)[, -c(1:2)]
      )

      # Write to CSV
      out_file <- paste0(
          path.3, "/meth.", region, ".", cntx, ".", treatment, ".csv"
      )
      write.csv(region_ctrl, out_file, row.names = FALSE)
  }

  
  cat("\n**\ndone :",treatment," - ",cntx,"\n**\n")
}

outer_results <- mclapply(c("mto1", "mto3"), function(treatment.loop) {
  inner_results <- mclapply(c("CG", "CHG", "CHH"), function(cntx.loop) {
    ave.meth.genes.levels(meth_GRlist[["wt"]], meth_GRlist[[treatment.loop]], paste0(treatment.loop, "_vs_wt"), cntx.loop)
  }, mc.cores = 3)
  return(inner_results)
}, mc.cores = 2)
#