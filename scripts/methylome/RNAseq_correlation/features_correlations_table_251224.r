library(dplyr)
library(GenomicRanges)
library(openxlsx)

RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag, DEG_log2FC = log2FoldChange) %>%
    dplyr::select(gene_id, DEG_log2FC, padj, pValue)

feature_file_fun <- function(context) {
    feature_file <- data.frame()
    for (ann in c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs", "TEG")) {

        # DMR results file
        ann_DMRs <- read.csv(paste0(
            "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/",
            context, "/", ann, "_", context, "_genom_annotations.csv"
        )) %>%
            select(
                gene_id, log2FC, context, type, Symbol, Computational_description
            ) %>%
            dplyr::rename(DMR_log2FC = log2FC)

        # correlation results file
        ann_corr <- read.csv(paste0(
            "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/Gene_feature/mto1/", context, "/", ann, ".corr.", context, ".mto1.csv"
        )) %>%
            select(-padj) %>%
            dplyr::rename(cor_pValue = pval)
        
        ann_merged = merge.data.frame(ann_corr, ann_DMRs, by = "gene_id", all.y = T)

        feature_file <- rbind(feature_file, ann_merged)
    }
    return(feature_file)
}

################
remove_dup_DMR <- function(y) {
    y <- as.character(unique(unlist(strsplit(y, ","))))
    paste(y, collapse = ",")
}
################

################
feture_res_table <- rbind(
    feature_file_fun("CG"),
    feature_file_fun("CHG"),
    feature_file_fun("CHH")
) %>% 
    mutate(
        DMR_log2FC = round(DMR_log2FC, 2),
        CG_DMRs = ifelse(grepl("CG", context), DMR_log2FC, NA),
        CHG_DMRs = ifelse(grepl("CHG", context), DMR_log2FC, NA),
        CHH_DMRs = ifelse(grepl("CHH", context), DMR_log2FC, NA)
    ) %>%
    merge.data.frame(RNAseq, ., by = "gene_id", all.y = TRUE) %>%
    mutate(tmp = paste(gene_id, type, sep = "_")) %>%
    group_by(tmp) %>%
    summarise(
        across(contains("CG_DMRs") | contains("CHG_DMRs") | contains("CHH_DMRs"), ~ remove_dup_DMR(paste(., collapse = ","))), # apply remove_dup_DMR function for DMR columns
        across(!contains("CG_DMRs") | contains("CHG_DMRs") | contains("CHH_DMRs"), dplyr::first) # for other columns
    ) %>%
    as.data.frame() %>%
    mutate(across(contains("_DMRs"), ~ gsub("NA", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",,", ",", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub("^,", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",$", "", .))) %>%
    mutate(across(contains("_DMRs"), ~ gsub(",", ", ", .))) %>%
    dplyr::relocate(CG_DMRs, CHG_DMRs, CHH_DMRs, .before = context) %>%
    dplyr::relocate(type, cor, cor_pValue, .after = gene_id) %>%
    select(-context, -DMR_log2FC, -tmp) %>%
    filter(pValue < 0.05) %>%
    filter(cor_pValue < 0.05) %>%
    arrange(pValue) %>%
    arrange(type) %>%
    mutate(across(contains("padj") | contains("pValue"), ~ gsub(" NA", NA, .)),
        pValue = as.numeric(formatC(.$pValue, format = "e", digits = 2)),
        padj = as.numeric(formatC(.$padj, format = "e", digits = 2)),
        DEG_log2FC = round(DEG_log2FC, 3),
        cor = round(cor, 3),
        cor_pValue = as.numeric(formatC(.$cor_pValue, format = "e", digits = 2))
    )


################
# for one sheet
xl_list <- list(feture_res_table = rbind(
    feture_res_table %>% filter(!type == "promoter"),
    feture_res_table %>% filter(type == "promoter")
))

################

################
clean_ASCII <- function(x) {
    x <- gsub("\002", " ", x)
    x <- gsub("\036", " ", x)
    return(x)
}
################
xl_headers <- names(xl_list[[1]])
DMRs_cols <- grep("_DMRs", xl_headers)
p_cols <- grep("cor_pValue|padj|^pValue$", xl_headers)
lfc_cols <- grep("^cor$|DEG_log2FC", xl_headers)
other_cols <- (grep("CHH_DMRs", xl_headers) + 2):length(xl_headers)
################
# save and edit EXCEL
wb <- createWorkbook()
# Define styles

# with border
if (F) {
    style_up <- createStyle(bgFill = "#f59d98", border = "TopBottomLeftRight", borderColour = "black")
    style_down <- createStyle(bgFill = "#c3ccf7", border = "TopBottomLeftRight", borderColour = "black")
    style_p <- createStyle(bgFill = "#f3e2c3", border = "TopBottomLeftRight", borderColour = "black")
    # style_other <- createStyle(bgFill = "#daf7d7", border = "TopBottomLeftRight", borderColour = "black")
    cell_n_font_style <- createStyle(border = "TopBottomLeftRight", borderColour = "black")
    header_style <- createStyle(textDecoration = "bold", border = "TopBottomLeftRight", borderStyle = "double")

    DMR_style_up <- createStyle(fgFill = "#f59d96", border = "TopBottomLeftRight", borderColour = "black")
    DMR_style_down <- createStyle(fgFill = "#c3ccf7", border = "TopBottomLeftRight", borderColour = "black")
    DMR_style_shared <- createStyle(fgFill = "#daf7d7", border = "TopBottomLeftRight", borderColour = "black")
}

# without border
if (T) {
    style_up <- createStyle(bgFill = "#f59d98", fontName = "Times New Roman")
    style_down <- createStyle(bgFill = "#c3ccf7", fontName = "Times New Roman")
    style_p <- createStyle(bgFill = "#f3e2c3", fontName = "Times New Roman")
    # style_other <- createStyle(bgFill = "#daf7d7", border = "TopBottomLeftRight", borderColour = "black")
    cell_n_font_style <- createStyle(fontName = "Times New Roman")
    header_style <- createStyle(textDecoration = "bold", border = "Bottom", borderStyle = "thick", fontName = "Times New Roman")

    DMR_style_up <- createStyle(fgFill = "#f59d96", fontName = "Times New Roman")
    DMR_style_down <- createStyle(fgFill = "#c3ccf7", fontName = "Times New Roman")
    DMR_style_shared <- createStyle(fgFill = "#daf7d7", fontName = "Times New Roman")
}

for (sheet_name in names(xl_list)) {
    df <- xl_list[[sheet_name]]
    df <- data.frame(lapply(df, clean_ASCII))
    df[, c(lfc_cols, p_cols)] <- sapply(df[, c(lfc_cols, p_cols)], as.numeric)

    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df)

    addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
    addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)


    # colors for DEGs values (minus as blue as plus as red)
    conditionalFormatting(wb, sheet_name, cols = lfc_cols[1], rows = 2:(nrow(df) + 1), rule = ">0", style = style_up)
    conditionalFormatting(wb, sheet_name, cols = lfc_cols[1], rows = 2:(nrow(df) + 1), rule = "<0", style = style_down)
    #
    conditionalFormatting(wb, sheet_name, cols = lfc_cols[2], rows = 2:(nrow(df) + 1), rule = ">0", style = style_up)
    conditionalFormatting(wb, sheet_name, cols = lfc_cols[2], rows = 2:(nrow(df) + 1), rule = "<0", style = style_down)

    # pValue columns
    conditionalFormatting(wb, sheet_name, cols = p_cols[1], rows = 2:(nrow(df) + 1), rule = "<0.05", style = style_p)
    conditionalFormatting(wb, sheet_name, cols = p_cols[2], rows = 2:(nrow(df) + 1), rule = "<0.05", style = style_p)
    conditionalFormatting(wb, sheet_name, cols = p_cols[3], rows = 2:(nrow(df) + 1), rule = "<0.05", style = style_p)

    # colors for DMRs values (minus as blue as plus as red)
    for (i.c in DMRs_cols) {
        for (i.r in 1:nrow(df)) {
            # for cells with both plus and minus direction
            if (length(grep("^-.*, [0-9]", df[i.r, i.c])) != 0 | length(grep("^[0-9].*, -[0-9]", df[i.r, i.c])) != 0) {
                addStyle(wb, sheet_name,
                    style = DMR_style_shared,
                    rows = i.r + 1, # first row in excel is the headers
                    cols = i.c,
                    gridExpand = FALSE
                )

                # for cells with minus direction
            } else if (length(grep("-", df[i.r, i.c])) != 0) {
                addStyle(wb, sheet_name,
                    style = DMR_style_down,
                    rows = i.r + 1, # first row in excel is the headers
                    cols = i.c,
                    gridExpand = FALSE
                )

                # for cells with plus direction
            } else if (length(grep("^[0-9]", df[i.r, i.c])) != 0) {
                addStyle(wb, sheet_name,
                    style = DMR_style_up,
                    rows = i.r + 1, # first row in excel is the headers
                    cols = i.c,
                    gridExpand = FALSE
                )
            }
        }
    }

    # the other columns
    # for (col in other_cols[other_cols %% 2 == 0]) {
    #    conditionalFormatting(wb, sheet_name, style = style_other, rule = "!=0", rows = 2:(nrow(df) + 1), cols = col, gridExpand = TRUE)
    #    conditionalFormatting(wb, sheet_name, style = style_other, rule = "==0", rows = 2:(nrow(df) + 1), cols = col, gridExpand = TRUE)
    # }

    # Remove gridlines
    showGridLines(wb, sheet_name, showGridLines = FALSE)
}
saveWorkbook(wb, paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/mto1_geneFeatures_correlation_table.xlsx"), overwrite = T)
