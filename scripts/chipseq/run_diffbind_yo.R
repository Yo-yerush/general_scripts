#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

sample_sheet <- args[1]
output_dir <- args[2]
test_condition <- args[3]
reference_condition <- args[4]
broad_peaks <- as.logical(args[5])

suppressPackageStartupMessages(library(DiffBind))

# Load samples and create a consensus peak set.
db <- dba(sampleSheet = sample_sheet)

if (broad_peaks) {
    db <- dba.count(db, minOverlap = 2, summits = FALSE)
} else {
    db <- dba.count(db, minOverlap = 2, summits = 250)
}

# Compare TEST against REFERENCE using DESeq2.
db <- dba.contrast(
    db,
    contrast = c("Condition", test_condition, reference_condition),
    minMembers = 2
)

db <- dba.analyze(
    db,
    method = DBA_DESEQ2,
    bBlacklist = FALSE,
    bGreylist = FALSE
)

saveRDS(db, file.path(output_dir, "diffbind_results.rds"))

# Export all differential regions.
all_regions <- dba.report(db, contrast = 1, method = DBA_DESEQ2, th = 1)
write.csv(
    as.data.frame(all_regions),
    file.path(output_dir, "differential_peaks_all.csv"),
    row.names = FALSE
)

# Export significant regions at FDR <= 0.05.
significant_regions <- dba.report(
    db, contrast = 1, method = DBA_DESEQ2, th = 0.05
)
write.csv(
    as.data.frame(significant_regions),
    file.path(output_dir, "differential_peaks_FDR_0.05.csv"),
    row.names = FALSE
)

# Diagnostic plots.
svg(file.path(output_dir, "PCA.svg"), width = 8, height = 7)
dba.plotPCA(db, attributes = DBA_CONDITION)
dev.off()

svg(file.path(output_dir, "correlation_heatmap.svg"), width = 8, height = 8)
plot(db)
dev.off()

svg(file.path(output_dir, "MA_plot.svg"), width = 8, height = 7)
dba.plotMA(db, contrast = 1, method = DBA_DESEQ2, th = 0.05)
dev.off()

svg(file.path(output_dir, "volcano_plot.svg"), width = 8, height = 7)
dba.plotVolcano(db, contrast = 1, method = DBA_DESEQ2, th = 0.05)
dev.off()

writeLines(
    c(
        paste("Contrast:", test_condition, "vs", reference_condition),
        paste("Positive Fold values indicate higher signal in", test_condition),
        paste("Significant regions:", length(significant_regions))
    ),
    file.path(output_dir, "summary.txt")
)

writeLines(
    capture.output(sessionInfo()),
    file.path(output_dir, "sessionInfo.txt")
)
