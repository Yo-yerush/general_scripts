library(dplyr)
library(writexl)
library(openxlsx)

###############
## fix it like mto1 paper
##############


des_file_1 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_description_file.csv.gz")

# des2 for protein.names and function columns
des_file_2 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/AT_description_file_230624.csv")
des_file_2 <- des_file_2[grep("^AT[0-5]G", des_file_2$locus_tag), c("locus_tag", "Protein.names", "Function..CC.")]
names(des_file_2)[c(1, 3)] <- c("gene_id", "Function")
des_file_2$Function <- gsub("FUNCTION:", "", des_file_2$Function)

des_file <- merge.data.frame(des_file_1, des_file_2, by = "gene_id", all = T) %>%
  relocate(Function, .before = note) %>%
  relocate(Protein.names, .after = Protein.families)

xl_list <- list()

for (treatment in c("mto1_vs_wt", "mto3_vs_wt", "dCGS_vs_EV", "SSE_high_vs_EV", "SSE_low_vs_EV", "SSE_high_vs_SSE_low")) {
  # ann.m = ifelse(ann == "genes", "Genes","Promoters")
  ############# methylome and RNAseq file filtered by corr files #############
  # corr_file = read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/corr_with_methylations/",treatment,"/",context,"/",ann,".corr.",context,".",treatment,".csv"))
  # corr_file = corr_file[corr_file$pval < 0.05, c("transcript_id","cor")]

  RNA_file <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/", treatment, "/all_genes_results_", treatment, ".csv"))


  RNA_file <- RNA_file[, c("gene_id", "log2FoldChange", "padj", "pValue")] # %>% filter(pValue < 0.05)
  names(RNA_file)[1:4] <- c("gene_id", "RNA_log2FC", "RNA_padj", "RNA_pvalue")
  RNA_file$RNA_log2FC <- round(RNA_file$RNA_log2FC, 3)
  RNA_file$RNA_padj <- as.numeric(sprintf("%.3e", RNA_file$RNA_padj))
  RNA_file$RNA_pvalue <- as.numeric(sprintf("%.3e", RNA_file$RNA_pvalue))

  # merge1 = merge.data.frame(corr_file,RNA_file, by = "transcript_id", all.y = T)
  meth_fun <- function(treatment.f, context.f, ann.f) {
    treatment.f <- ifelse(grepl("SSE_H_vs_EV", treatment.f), "SSE_high_vs_EV", treatment.f)
    treatment.f <- ifelse(grepl("SSE_L_vs_EV", treatment.f), "SSE_low_vs_EV", treatment.f)

    meth_file.f <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/", treatment.f, "/genome_annotation/", context.f, "/", ann.f, "_", context.f, "_genom_annotations.csv"))

    # meth_file.f = read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment.f,"_vs_wt/genome_annotation/",context.f,"/",ann.f,"_",context.f,"_genom_annotations.csv"))
    meth_file.f$log2FoldChange <- round(log2(meth_file.f$proportion2 / meth_file.f$proportion1), 3)
    # meth_file.f = meth_file.f[,c("gene_id","seqnames","start","end","log2FoldChange")]
    meth_file.f <- meth_file.f[, c("gene_id", "log2FoldChange")]
    names(meth_file.f)[2] <- paste0(context.f, "_", ann.f)
    return(meth_file.f)
  }


  meth_list <- list()
  if (treatment == "SSE_high_vs_SSE_low") {
    meth_list[["CG_genes"]] <- data.frame(gene_id = NA, CG_Genes = NA)
    meth_list[["CG_promoters"]] <-  data.frame(gene_id = NA, CG_Promoters = NA)
  } else {
    meth_list[["CG_genes"]] <- meth_fun(treatment, "CG", "Genes")
    meth_list[["CG_promoters"]] <- meth_fun(treatment, "CG", "Promoters")
  }
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
  merge3 <- merge3[order(merge3$RNA_padj), ]
  # merge3 <- merge3[order(merge3$cor, decreasing = T), ]

  xl_list[[treatment]] <- merge3
}


########################################################

write_xlsx(xl_list, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_genes.xlsx")

########################################################
# its not working, maybe beacouse of a file zise issue
if (F) {
  # save and edit EXCEL
  wb <- createWorkbook()
  # Define styles
  style_up <- createStyle(bgFill = "#f59d98")
  style_down <- createStyle(bgFill = "#c3ccf7")
  style_p <- createStyle(bgFill = "#f7deb0")
  style_other <- createStyle(bgFill = "#daf7d7")
  cell_n_font_style <- createStyle(border = "TopBottomLeftRight", borderColour = "black")
  header_style <- createStyle(textDecoration = "bold", border = "TopBottomLeftRight", borderStyle = "double")

  for (sheet_name in names(xl_list)) {
    addWorksheet(wb, sheet_name)
    df <- xl_list[[sheet_name]]
    writeData(wb, sheet_name, df)

    conditionalFormatting(wb, sheet_name, cols = grep("RNA_padj", names(df)), rows = 2:(nrow(df) + 1), rule = "<0.05", style = style_p)

    quantitative_cols <- grep("RNA_log2FC|CG_|CHG_|CHH_", names(df))
    for (col in quantitative_cols) {
      conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df) + 1), rule = ">0", style = style_up)
      conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df) + 1), rule = "<0", style = style_down)
    }

    other_cols <- (max(quantitative_cols) + 1):ncol(df)
    for (col in other_cols[other_cols %% 2 == 1]) {
      conditionalFormatting(wb, sheet_name, style = style_other, rule = "!=0", rows = 2:(nrow(df) + 1), cols = col, gridExpand = TRUE)
      conditionalFormatting(wb, sheet_name, style = style_other, rule = "==0", rows = 2:(nrow(df) + 1), cols = col, gridExpand = TRUE)
    }


    addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
    addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  }
  saveWorkbook(wb, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_transcripts_040324.xlsx", overwrite = T)
}

#
