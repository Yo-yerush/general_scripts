windowSize <- function(x, mcol_name, windowSize = 1.5e5) {

    # chromosome lengths
    seqlens <- vapply(split(end(x), seqnames(x)), max, numeric(1))
    seqlengths(x) <- seqlens[seqlevels(x)]

    # windows by window size
    windows <- tileGenome(seqlens,
        tilewidth = windowSize,
        cut.last.tile.in.chrom = TRUE
    )

    # map to windows and calculate mean value
    hits <- findOverlaps(windows, x, ignore.strand = TRUE)
    mValue <- tapply(mcols(x)[[mcol_name]][subjectHits(hits)],
        queryHits(hits),
        mean,
        na.rm = TRUE
    )
    mcols(windows)$mean_value <- NA_real_
    mcols(windows)$mean_value[as.integer(names(mValue))] <- mValue

    windows
}
