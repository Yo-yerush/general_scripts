# YO - 100725
# average methylation within each gene features
# this output can use for transcriptome-methylome correlation analysis
# will run **42** jobs in parallel (max)

average_methylation_over_gene_features <- function(sample_file_path, output_path = "./") {
    library(DMRcaller)
    library(parallel)
    library(dplyr)

    source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/load_replicates.R")
    source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/trimm_and_rename_seq.R")

    # read annotation file from github
    annotations_url <- "https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/annotation_files/Methylome.At_annotations.csv.gz"
    annotation_file <- read.csv(
        gzcon(url(annotations_url, open = "rb"), text = TRUE),
        encoding = "UTF-8"
    )

    ##########################################
    # define varibals
    sample_file <- read.table(sample_file_path, sep = "\t")
    var1_name <- unique(sample_file[, 1])[1]
    var2_name <- unique(sample_file[, 1])[2]
    var1_path <- sample_file[sample_file[, 1] == var1_name, 2]
    var2_path <- sample_file[sample_file[, 1] == var2_name, 2]

    var_args <- list(
        list(path = var1_path, name = var1_name),
        list(path = var2_path, name = var2_name)
    )

    # load
    n_cores_in_fun <- max(c(length(var1_path), length(var2_path)))
    load_vars <- mclapply(var_args, function(x) load_replicates(x$path, n_cores_in_fun, x$name),
        mc.cores = 2
    )

    load_var1 <- trimm_and_rename(load_vars[[1]]$methylationDataReplicates)
    load_var2 <- trimm_and_rename(load_vars[[2]]$methylationDataReplicates)

    # add ratio column for each sample (methylated / unmethylated)
    ratio_cols <- function(x, x_name) {
        n_samples <- length(grep("readsM", names(mcols(x))))
        x <- as.data.frame(x)
        for (i in 1:n_samples) {
            ratio_col_name <- paste0(x_name, ".", i)
            x[ratio_col_name] <- x[paste0("readsM", i)] / x[paste0("readsN", i)]
        }
        x <- makeGRangesFromDataFrame(x, keep.extra.columns = T)
        return(x)
    }
    meth_var1 <- ratio_cols(load_var1, var1_name)
    meth_var2 <- ratio_cols(load_var2, var2_name)

    meth_GRlist <- list(meth_var1, meth_var2)
    names(meth_GRlist) <- c(var1_name, var2_name)

    ##########################################

    ave.meth.genes.levels <- function(meth_ctrl, meth_trmt, ctrl_name, trmt_name, cntx) {
        treatment <- paste0(trmt_name, "_vs_", ctrl_name)
        path.0 <- output_path
        path.1 <- paste0(path.0, "average.meth.genes.levels")
        path.2 <- paste0(path.1, "/", treatment)
        path.3 <- paste0(path.2, "/", cntx)
        dir.create(path.1, showWarnings = F)
        dir.create(path.2, showWarnings = F)
        dir.create(path.3, showWarnings = F)

        meth_ctrl_cntx <- meth_ctrl[which(meth_ctrl$context == cntx)]
        meth_trmt_cntx <- meth_trmt[which(meth_trmt$context == cntx)]

        meth_ctrl_cntx <- meth_ctrl_cntx[, grep(ctrl_name, names(meth_ctrl_cntx@elementMetadata))]
        meth_trmt_cntx <- meth_trmt_cntx[, grep(strsplit(trmt_name, "_")[[1]], names(meth_trmt_cntx@elementMetadata))]

        # Specify the annotation types
        annotation_types <- c("promoter", "gene", "CDS", "intron", "five_prime_UTR", "three_prime_UTR", "transposable_element_gene")

        # Read in each annotation type in a loop
        dmrs_list <- lapply(annotation_types, function(region) {
            if (region == "promoter") {
                annotation_file %>%
                    filter(gene_model_type == "protein_coding") %>%
                    filter(type == "gene") %>%
                    dplyr::select(-gene_model_type) %>%
                    makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>%
                    promoters(., upstream = 2000, downstream = 0)
                # } else if (region == "transposable_element_gene") {
                #    annotation_file %>%
                #        filter(type == "transposable_element_gene") %>%
                #        dplyr::select(-gene_model_type) %>%
                #        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
            } else {
                annotation_file %>%
                    filter(type == region) %>%
                    dplyr::select(-gene_model_type) %>%
                    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
            }
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
        meth.list <- mclapply(seq_along(dmrs_list), function(DL) {
            mclapply(seq_along(meth_list),
                function(ML) fun4parallel(ML, dmrs_list, DL),
                mc.cores = 2
            )
        }, mc.cores = length(dmrs_list))

        # final loop:
        for (i in seq_along(annotation_types)) {
            region <- annotation_types[i]

            # Control & treatment are each a GRanges object
            region_ctrl <- meth.list[[i]][[1]]
            region_trmt <- meth.list[[i]][[2]]

            # Combine control & treatment columns (minus the first two columns, as in your original code)
            mcols(region_ctrl) <- cbind(
                mcols(region_ctrl),
                mcols(region_trmt)[, -grep("gene_id|type", names(mcols(region_trmt)))]
            )

            # Write to CSV
            out_file <- paste0(
                path.3, "/meth.", region, ".", cntx, ".", treatment, ".csv"
            )
            write.csv(region_ctrl, out_file, row.names = FALSE)
        }


        cat("\n**\ndone :", treatment, " - ", cntx, "\n**\n")
    }


    inner_results <- mclapply(c("CG", "CHG", "CHH"), function(cntx.loop) {
        ave.meth.genes.levels(meth_GRlist[[var1_name]], meth_GRlist[[var2_name]], var1_name, var2_name, cntx.loop)
    }, mc.cores = 3)
}
