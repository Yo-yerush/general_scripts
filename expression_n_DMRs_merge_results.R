expression_n_DMRs_merge_results <- function(
    treatment,
    DMRs_results_directory = "/PATH/TO/Methylome.At/results/mut_vs_wt/genome_annotation/",
    RNAseq_results_directory = "/PATH/TO/rnaseq_results/mut_vs_wt/") {

  des_file <- data.table::fread("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/refs/heads/main/annotation_files/Methylome.At_description_file.csv.gz")

  RNA_file <- read.csv(paste0(RNAseq_results_directory, "all_genes_results_", treatment, ".csv"))
  RNA_file <- RNA_file[, c("gene_id", "log2FoldChange", "padj", "pValue")] # %>% filter(pValue < 0.05)
  names(RNA_file)[1:4] <- c("gene_id", "RNA_log2FC", "RNA_padj", "RNA_pvalue")
  RNA_file$RNA_log2FC <- round(RNA_file$RNA_log2FC, 3)
  RNA_file$RNA_padj <- as.numeric(sprintf("%.3e", RNA_file$RNA_padj))
  RNA_file$RNA_pvalue <- as.numeric(sprintf("%.3e", RNA_file$RNA_pvalue))

  # merge1 = merge.data.frame(corr_file,RNA_file, by = "transcript_id", all.y = T)
  meth_fun <- function(treatment.f, context.f, ann.f) {
    tryCatch(
      {
        meth_file.f <- read.csv(paste0(DMRs_results_directory, context.f, "/", ann.f, "_", context.f, "_genom_annotations.csv"))
        meth_file.f$log2FoldChange <- round(log2(meth_file.f$proportion2 / meth_file.f$proportion1), 3)
        meth_file.f <- meth_file.f[, c("gene_id", "log2FoldChange")]
      },
      error = function(e) {
        meth_file.f <- data.frame(gene_id = NA, log2FoldChange = NA)
      }
    )
    names(meth_file.f)[2] <- paste0(context.f, "_", ann.f)
    return(meth_file.f)
  }

  meth_list <- list()
  meth_list[["CG_genes"]] <- meth_fun(treatment, "CG", "Genes")
  meth_list[["CG_promoters"]] <- meth_fun(treatment, "CG", "Promoters")
  meth_list[["CHG_genes"]] <- meth_fun(treatment, "CHG", "Genes")
  meth_list[["CHH_genes"]] <- meth_fun(treatment, "CHH", "Genes")
  meth_list[["CHG_promoters"]] <- meth_fun(treatment, "CHG", "Promoters")
  meth_list[["CHH_promoters"]] <- meth_fun(treatment, "CHH", "Promoters")

  merged_meth <- meth_list[[1]]
  for (i in 2:length(meth_list)) {
    merged_meth <- merge.data.frame(merged_meth, meth_list[[i]], by = "gene_id", all = TRUE)
  }

  merge2 <- merge.data.frame(RNA_file, merged_meth, by = "gene_id", all.x = T)
  merge3 <- merge.data.frame(merge2, des_file, by = "gene_id", all.x = T)

  merge3[order(merge3$RNA_padj), ]
}
