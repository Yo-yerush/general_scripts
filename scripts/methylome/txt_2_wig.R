# For files with format: Chr1 - 1 (chromosome, strand, position)

txt_to_wig <- function(input_file, output_file = "output.bw", smRNA = FALSE, h3_file = NULL, bigwig = TRUE, bin_size = 50) {

    # Check if output file exists
    if (file.exists(output_file)) {
        cat("Output file already exists:", output_file, "\n")
        return(invisible())
    }
    
    # disable scientific notation
    options(scipen = 999)

    tair10_sizes <- c(
        Chr1 = 30427671, Chr2 = 19698289, Chr3 = 23459830,
        Chr4 = 18585056, Chr5 = 26975502, ChrC = 154478, ChrM = 367808 # 366924
    )

    cat("Reading", basename(input_file), "\n")

    # if (smRNA) {
    #     seq_df <- read.table(input_file, sep = "\t", stringsAsFactors = FALSE)
    #     colnames(seq_df) <- c("chr", "strand", "position", "cov")
    #     seq_df$cov <- as.numeric(seq_df$cov)
    # } else {
    #     seq_df <- read.table(input_file, sep = " ", stringsAsFactors = FALSE)
    #     colnames(seq_df) <- c("chr", "strand", "position")
    # }

    seq_df <- read.table(input_file, sep = ifelse(smRNA, "\t", " "), stringsAsFactors = FALSE)
    colnames(seq_df) <- c("chr", "strand", "position")

    cat("Processing", nrow(seq_df), "reads\n")

    seq_df$position <- as.numeric(seq_df$position)
    # seq_df <- seq_df[!is.na(seq_df$position) & seq_df$position >= 0, ]

    if (!is.null(h3_file)) {
        H3 <- read.table(h3_file, sep = " ", stringsAsFactors = FALSE)
        colnames(H3) <- c("chr", "strand", "position")

        # keep unmatched rows from H3Kn file and H3
        seq_df <- seq_df[!paste(seq_df$chr, seq_df$position) %in% paste(H3$chr, H3$position), ]

        cat("H3 filtered:", nrow(seq_df), "reads\n")
    }

    # Count reads per position for each chromosome
    coverage_list <- list()

    for (chr in unique(seq_df$chr)) {
        chr_data <- seq_df[seq_df$chr == chr, ]
        pos_counts <- table(chr_data$position)
        # if (smRNA) {
        #     cov_counts <- chr_data$cov
        # }

        coverage_list[[chr]] <- data.frame(
            position = as.numeric(names(pos_counts)),
            # count = ifelse(smRNA, as.numeric(cov_counts), as.numeric(pos_counts))
            count = ifelse(smRNA, as.numeric(pos_counts), as.numeric(pos_counts))
        )
    }
    #####################################
    # Write WIG file
    if (!bigwig) {
        cat("Writing", output_file, "\n")
        con <- file(output_file, "w")
        writeLines("track type=wiggle_0", con)
        for (chr in sort(names(coverage_list))) {
            if (nrow(coverage_list[[chr]]) > 0) {
                writeLines(paste0("variableStep chrom=", chr), con)
                chr_cov <- coverage_list[[chr]][order(coverage_list[[chr]]$position), ]
                for (i in 1:nrow(chr_cov)) {
                    writeLines(paste(chr_cov$position[i], chr_cov$count[i]), con)
                }
            }
        }
        close(con)
        cat("WIG file created:", output_file, "\n")
    } else {
        # ---- BigWig branch (direct write) ----
        # Build GRanges of 1-bp bins with a 'score' column
        gr_list <- lapply(names(coverage_list), function(i_chr) {
            df <- coverage_list[[i_chr]]
            if (!nrow(df)) {
                return(GenomicRanges::GRanges())
            }
            df <- df[order(df$position), , drop = FALSE]
            GenomicRanges::GRanges(
                seqnames = i_chr,
                ranges = IRanges::IRanges(start = df$position, width = 1),
                score = df$count,
                seqinfo = GenomeInfoDb::Seqinfo(
                    seqnames   = names(tair10_sizes),
                    seqlengths = as.numeric(tair10_sizes),
                    genome     = "TAIR10"
                )
            )
        })
        gr <- do.call(c, gr_list)

        # # Keep only canonical TAIR10 seqlevels and attach lengths
        # gr <- GenomicRanges::keepSeqlevels(
        #     gr,
        #     intersect(GenomeInfoDb::seqlevels(gr), names(tair10_sizes)),
        #     pruning.mode = "coarse"
        # )
        GenomeInfoDb::seqlengths(gr) <- tair10_sizes[GenomeInfoDb::seqlevels(gr)]
        GenomeInfoDb::genome(gr) <- "TAIR10"

        # compute binned density
        cov_rle <- GenomicRanges::coverage(gr, weight = "score")
        bins <- GenomicRanges::tileGenome(
            seqlengths = tair10_sizes,
            tilewidth = bin_size,
            cut.last.tile.in.chrom = TRUE
        )
        binned <- GenomicRanges::binnedAverage(bins, cov_rle, "score")


        # Derive .bw name from output_file if needed
        bw_file <- if (grepl("\\.bw$", output_file, ignore.case = TRUE)) {
            output_file
        } else {
            sub("\\.wig(\\.gz)?$", ".bw", output_file, ignore.case = TRUE)
        }
        if (!grepl("\\.bw$", bw_file, ignore.case = TRUE)) bw_file <- paste0(bw_file, ".bw")

        cat("Writing BigWig:", bw_file, "\n")
        rtracklayer::export(binned, bw_file, format = "BigWig")
        cat("BigWig created:", bw_file, "\n")
    }
}

# input_path <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/GSE51304_RAW/GSM1667165_WT_H3K9me2_ChIP_rep_processed.txt.gz"
# input_h3_path <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/GSE51304_RAW/GSM1667164_WT_H3_ChIP_rep_processed.txt.gz"
# output_path <- "WT_H3K9me2_ChIP_rep_processed.wig.gz"
# txt_to_wig(input_path, output_path, smRNA = F, h3_file = input_h3_path)
