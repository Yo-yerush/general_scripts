windowSize <- function(x, mcol_name, windowSize = 1.5e5, mean_value = T, count_sites = T, sum_value = F) {
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


    if (mean_value) {
        mValue <- tapply(mcols(x)[[mcol_name]][subjectHits(hits)],
            queryHits(hits),
            mean,
            na.rm = TRUE
        )
        mcols(windows)$mean_value <- NA_real_
        mcols(windows)$mean_value[as.integer(names(mValue))] <- mValue
    }

    if (count_sites) {
        windows$n_sites <- tabulate(queryHits(hits), nbins = length(windows))
    }

    if (sum_value) {
        sValue <- tapply(mcols(x)[[mcol_name]][subjectHits(hits)],
            queryHits(hits),
            sum,
            na.rm = TRUE
        )
        mcols(windows)$sum_value <- NA_real_
        mcols(windows)$sum_value[as.integer(names(sValue))] <- sValue
    }

    windows
}
