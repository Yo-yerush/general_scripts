expression_n_DMRs_merge_results <- function(
    DMRs_results_directory,
    RNAseq_results,
    designed_excel_output = NULL, # if not NULL, write the output file prefix
    make_it_unique = FALSE, # eill make unique IDs data frame but with non-numeric values got the DMRs log2FC
    DMR_sep = ", ", # if unique
    empty_DMR = "" # if unique
    ) {
  des_file <- data.table::fread("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/refs/heads/main/annotation_files/Methylome.At_description_file.csv.gz")

  if (any(duplicated(RNA_file[, 1]))) stop("duplicated IDs (column 1)")
  if (!all(apply(RNA_file[, 2:4], 2, is.numeric))) stop("non-numeric values (columns 2-4)")

  RNA_file <- RNAseq_results[, 1:4]
  RNA_new_names <- c("gene_id", "RNA_log2FC", "RNA_padj", "RNA_pvalue")
  names(RNA_file) <- RNA_new_names
  RNA_file$RNA_log2FC <- round(RNA_file$RNA_log2FC, 3)
  RNA_file$RNA_padj <- as.numeric(sprintf("%.3e", RNA_file$RNA_padj))
  RNA_file$RNA_pvalue <- as.numeric(sprintf("%.3e", RNA_file$RNA_pvalue))

  meth_fun <- function(context.f, ann.f) {
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
  meth_list[["CG_genes"]] <- meth_fun("CG", "Genes")
  meth_list[["CG_promoters"]] <- meth_fun("CG", "Promoters")
  meth_list[["CHG_genes"]] <- meth_fun("CHG", "Genes")
  meth_list[["CHH_genes"]] <- meth_fun("CHH", "Genes")
  meth_list[["CHG_promoters"]] <- meth_fun("CHG", "Promoters")
  meth_list[["CHH_promoters"]] <- meth_fun("CHH", "Promoters")

  merged_meth <- meth_list[[1]]
  for (i in 2:length(meth_list)) {
    merged_meth <- merge.data.frame(merged_meth, meth_list[[i]], by = "gene_id", all = TRUE)
  }

  merged_all <- merge.data.frame(RNA_file, merged_meth, by = "gene_id", all.x = T)
  merged_out <- merge.data.frame(merged_all, des_file, by = "gene_id", all.x = T)
  merged_out <- merged_out[order(merged_out$RNA_padj), ]

  if (make_it_unique) {
    merged_out <- unique_DMRs_columns(merged_out, DMR_sep = DMR_sep, empty_DMR = empty_DMR)
  }

  if (!is.null(designed_excel_output)) {
    openxlsx::saveWorkbook(
      merged_excel_formate(merged_out, designed_excel_output),
      paste0(designed_excel_output, ".xlsx"),
      overwrite = T
    )
  }

  return(merged_out)
}

##############################################################

remove_dup_DMR <- function(y) {
  y <- as.character(unique(unlist(strsplit(y, ","))))
  paste(y, collapse = ",")
}

##############################################################

clean_ASCII <- function(x) {
  x <- gsub("\001", " ", x)
  x <- gsub("\002", " ", x)
  x <- gsub("\036", " ", x)
  return(x)
}

##############################################################

unique_DMRs_columns <- function(x, DMR_sep = ", ", empty_DMR = "") {
  suppressMessages({
    library(dplyr)
    library(stringr)
  })

  x %>%
    group_by(gene_id) %>%
    summarise(
      across(contains("_Genes") | contains("_Promoters"), ~ remove_dup_DMR(paste(., collapse = ","))), # apply remove_dup_DMR function for DMR columns
      across(!contains("_Genes") & !contains("_Promoters"), first) # for other columns
    ) %>%
    as.data.frame() %>%
    mutate(across(contains("_Genes") | contains("_Promoters"), ~ str_replace_all(.x, "NA", empty_DMR))) %>%
    mutate(across(contains("_Genes") | contains("_Promoters"), ~ str_replace_all(.x, ",", DMR_sep))) %>% # change delimiter
    relocate(contains("_Genes") | contains("_Promoters"), .after = RNA_pvalue) %>%
    arrange(RNA_padj)
}

##############################################################

merged_excel_formate <- function(df, sheet_name) {
  library(openxlsx)

  numeric_cols <- grep("RNA_", names(df))
  p_cols <- grep("RNA_p", names(df))
  lfc_rna_cols <- grep("RNA_log2FC", names(df))
  lfc_dmr_cols <- grep("CG_|CHG_|CHH_", names(df))
  other_cols <- (grep("CHH_Promoters", names(df)) + 1):length(names(df))
  ################
  # save and edit EXCEL
  wb <- createWorkbook()
  # Define styles
  style_up <- createStyle(fontName = "Times New Roman", bgFill = "#f59d98", fgFill = "#f59d98", border = "TopBottomLeftRight", borderColour = "black")
  style_down <- createStyle(fontName = "Times New Roman", bgFill = "#c3ccf7", fgFill = "#c3ccf7", border = "TopBottomLeftRight", borderColour = "black")
  style_p <- createStyle(fontName = "Times New Roman", bgFill = "#f7deb0", border = "TopBottomLeftRight", borderColour = "black")
  style_other <- createStyle(fontName = "Times New Roman", bgFill = "#daf7d7", fgFill = "#daf7d7", border = "TopBottomLeftRight", borderColour = "black")
  cell_n_font_style <- createStyle(fontName = "Times New Roman", border = "TopBottomLeftRight", borderColour = "black")
  header_style <- createStyle(fontName = "Times New Roman", textDecoration = "bold", border = "TopBottomLeftRight", borderStyle = "double")

  addWorksheet(wb, sheet_name)

  df <- data.frame(lapply(df, clean_ASCII))
  df[, numeric_cols] <- sapply(df[, numeric_cols], as.numeric)

  writeData(wb, sheet_name, df)

  addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
  addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)

  conditionalFormatting(wb, sheet_name, cols = p_cols[1], rows = 2:(nrow(df) + 1), rule = "<0.05", style = style_p)
  conditionalFormatting(wb, sheet_name, cols = p_cols[2], rows = 2:(nrow(df) + 1), rule = "<0.05", style = style_p)

  for (col in lfc_rna_cols) {
    conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df) + 1), rule = ">0", style = style_up)
    conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df) + 1), rule = "<0", style = style_down)
  }

  # colors for DMRs values (minus as blue as plus as red)
  for (i.c in lfc_dmr_cols) {
    for (i.r in 1:nrow(df)) { # first row in excel is the headers
      # for cells with both plus and minus direction
      if (length(grep("^-.*, [0-9]", df[i.r, i.c])) != 0 | length(grep("^[0-9].*, -[0-9]", df[i.r, i.c])) != 0) {
        addStyle(wb, sheet_name,
          style = style_other,
          rows = i.r + 1,
          cols = i.c,
          gridExpand = TRUE
        )

        # for cells with minus direction
      } else if (length(grep("-", df[i.r, i.c])) != 0) {
        addStyle(wb, sheet_name,
          style = style_down,
          rows = i.r + 1,
          cols = i.c,
          gridExpand = TRUE
        )

        # for cells with plus direction
      } else if (length(grep("^[0-9]", df[i.r, i.c])) != 0) {
        addStyle(wb, sheet_name,
          style = style_up,
          rows = i.r + 1,
          cols = i.c,
          gridExpand = TRUE
        )
      }
    }
  }

  for (col in other_cols[other_cols %% 2 == 0]) {
    conditionalFormatting(wb, sheet_name, style = style_other, rule = "!=0", rows = 2:(nrow(df) + 1), cols = col, gridExpand = TRUE)
    conditionalFormatting(wb, sheet_name, style = style_other, rule = "==0", rows = 2:(nrow(df) + 1), cols = col, gridExpand = TRUE)
  }

  return(wb)
}
