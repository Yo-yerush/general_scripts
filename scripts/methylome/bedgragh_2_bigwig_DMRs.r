# Load required library
library(rtracklayer)
library(org.At.tair.db)

output_dir <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/bigWig_files_results/"
results_dir <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/"

DMRs_2_bigWig <- function(var1, var2, context) {
    DMRsReplicates <- read.csv(paste0(results_dir, var2, "_vs_", var1, "/DMRs_", context, "_", var2, "_vs_", var1, ".csv"))
    file_name <- paste(var2, context, sep = "_")

    # keep genome location and methylation proportion
    DMRsReplicates_wig <- DMRsReplicates[, c("seqnames", "start", "end")]
    DMRsReplicates_wig$score <- round(log2(DMRsReplicates$proportion2 / DMRsReplicates$proportion1), 2)
    DMRsReplicates_wig <- makeGRangesFromDataFrame(DMRsReplicates_wig, keep.extra.columns = T)

    # get seqlength for arabidopsis
    chrInfo <- as.list(org.At.tairCHRLENGTHS)
    seqlengths <- unlist(chrInfo)[1:5]
    names(seqlengths) <- paste0("Chr", 1:5)

    # save as BigWig file
    seqlengths(DMRsReplicates_wig) <- seqlengths
    export(DMRsReplicates_wig, paste0(output_dir, file_name, ".bw"), format = "bigWig")
}

dir.create(output_dir, showWarnings = F)

for (cntx in c("CG", "CHG", "CHH")) {
    DMRs_2_bigWig("wt", "mto1", cntx)
    DMRs_2_bigWig("wt", "mto3", cntx)
    DMRs_2_bigWig("EV", "dCGS", cntx)
    DMRs_2_bigWig("EV", "SSE_high", cntx)
    DMRs_2_bigWig("EV", "SSE_low", cntx)
}
