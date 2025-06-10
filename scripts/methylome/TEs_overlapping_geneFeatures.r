library(dplyr)
library(GenomicRanges)
library(openxlsx)

RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
    dplyr::rename(gene_id = locus_tag, DEG_log2FC = log2FoldChange) %>%
    dplyr::select(gene_id, DEG_log2FC, padj, pValue)

overlapped_TE_context <- function(context) {
    TE <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Transposable_Elements_", context, "_genom_annotations.csv")) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    overlapped <- data.frame()
    for (i.feature in c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs")) {
        gene_feature <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/", i.feature, "_", context, "_genom_annotations.csv")) %>%
            makeGRangesFromDataFrame(keep.extra.columns = TRUE)

        m <- findOverlaps(gene_feature, TE)
        m_df <- gene_feature[queryHits(m)]
        m_df$overlapping_TE <- paste(TE$Transposon_Super_Family[subjectHits(m)], TE$Transposon_Family[subjectHits(m)], sep = ", ")
        #m_df$Transposon_Super_Family <- TE$Transposon_Super_Family[subjectHits(m)]
        #m_df$Transposon_Family <- TE$Transposon_Family[subjectHits(m)]

        overlapped <- rbind(overlapped, as.data.frame(m_df))
    }

    return(overlapped)
}

################
remove_dup_DMR <- function(y) {
    y <- as.character(unique(unlist(strsplit(y, ","))))
    paste(y, collapse = ",")
}
################

################
overlapped_TE <- rbind(
    overlapped_TE_context("CG"),
    overlapped_TE_context("CHG"),
    overlapped_TE_context("CHH")
) %>%
    select(
        -seqnames, -start, -end, -width, -strand, -pValue, -direction, -regionType,
        -sumReadsM1, -sumReadsN1, -proportion1, -sumReadsM2, -sumReadsN2, -proportion2, -cytosinesCount
    ) %>%
    dplyr::rename(DMR_log2FC = log2FC) %>%
    mutate(
        DMR_log2FC = round(DMR_log2FC, 3),
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
    dplyr::relocate(type, overlapping_TE, .after = gene_id) %>%
    dplyr::relocate(Computational_description, .before = Short_description) %>%
    select(-context, -DMR_log2FC, -tmp, -(Short_description : last_col())) %>%
    filter(pValue < 0.05) %>%
    arrange(pValue) %>%
    arrange(type) %>%
    mutate(across(contains("padj") | contains("pValue"), ~ gsub(" NA", NA, .)),
        pValue = as.numeric(formatC(.$pValue, format = "e", digits = 3)),
        padj = as.numeric(formatC(.$padj, format = "e", digits = 3)),
        DEG_log2FC = round(DEG_log2FC, 3)
    )


################
# for one sheet
xl_list <- list(overlapped_TE = rbind(
    overlapped_TE %>% filter(!type == "promoter"),
    overlapped_TE %>% filter(type == "promoter")
))

# Create a dataframe with counts of upregulated and downregulated DEGs and DMRs - FOR RACHEL!!!!
gene_feature_counts <- xl_list[[1]] %>%
    mutate(
        across(contains("_DMRs"), ~ as.numeric(gsub(",.*", "", .)))
    ) %>%
    filter(pValue < 0.05) %>%
    mutate(
        upregulated_DEGs = ifelse(DEG_log2FC > 0, 1, 0),
        downregulated_DEGs = ifelse(DEG_log2FC < 0, 1, 0),
        CG_hyper_DMRs = ifelse(CG_DMRs > 0, 1, 0),
        CG_hypo_DMRs = ifelse(CG_DMRs < 0, 1, 0),
        CHG_hyper_DMRs = ifelse(CHG_DMRs > 0, 1, 0),
        CHG_hypo_DMRs = ifelse(CHG_DMRs < 0, 1, 0),
        CHH_hyper_DMRs = ifelse(CHH_DMRs > 0, 1, 0),
        CHH_hypo_DMRs = ifelse(CHH_DMRs < 0, 1, 0)
    ) %>%
    group_by(type) %>%
    summarise(
        upregulated_DEGs = sum(upregulated_DEGs, na.rm = TRUE),
        downregulated_DEGs = sum(downregulated_DEGs, na.rm = TRUE),
        CG_hyper_DMRs = sum(CG_hyper_DMRs, na.rm = TRUE),
        CG_hypo_DMRs = sum(CG_hypo_DMRs, na.rm = TRUE),
        CHG_hyper_DMRs = sum(CHG_hyper_DMRs, na.rm = TRUE),
        CHG_hypo_DMRs = sum(CHG_hypo_DMRs, na.rm = TRUE),
        CHH_hyper_DMRs = sum(CHH_hyper_DMRs, na.rm = TRUE),
        CHH_hypo_DMRs = sum(CHH_hypo_DMRs, na.rm = TRUE)
    ) %>%
    as.data.frame() %>%
    write.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_overlapping_geneFeatures/mto1_direction_count_TEs_overlapping_geneFeatures.csv", row.names = F)

# Add the new dataframe to the list of sheets
xl_list$gene_feature_counts <- gene_feature_counts

# for two sheets
#xl_list <- list(
#    Gene_bodies = overlapped_TE %>% filter(!type == "promoter"),
#    Promoters = overlapped_TE %>% filter(type == "promoter")
#)
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
p_cols <- grep("padj|pValue", xl_headers)
lfc_cols <- grep("DEG_log2FC", xl_headers)
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
    style_up <- createStyle(bgFill = "#f59d98")
    style_down <- createStyle(bgFill = "#c3ccf7")
    style_p <- createStyle(bgFill = "#f3e2c3")
    # style_other <- createStyle(bgFill = "#daf7d7", border = "TopBottomLeftRight", borderColour = "black")
    cell_n_font_style <- createStyle(fontName = "Times New Roman")
    header_style <- createStyle(textDecoration = "bold", border = "Bottom", borderStyle = "thick")

    DMR_style_up <- createStyle(fgFill = "#f59d96")
    DMR_style_down <- createStyle(fgFill = "#c3ccf7")
    DMR_style_shared <- createStyle(fgFill = "#daf7d7")
}

for (sheet_name in names(xl_list)) {
    df <- xl_list[[sheet_name]]
    df <- data.frame(lapply(df, clean_ASCII))
    df[, c(lfc_cols, p_cols)] = sapply(df[, c(lfc_cols, p_cols)], as.numeric)

    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df)

    addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
    addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)


    # colors for DEGs values (minus as blue as plus as red)
    conditionalFormatting(wb, sheet_name, cols = lfc_cols, rows = 2:(nrow(df) + 1), rule = ">0", style = style_up)
    conditionalFormatting(wb, sheet_name, cols = lfc_cols, rows = 2:(nrow(df) + 1), rule = "<0", style = style_down)

    # pValue columns
     conditionalFormatting(wb, sheet_name, cols = p_cols[1], rows = 2:(nrow(df) + 1), rule = "<0.05", style = style_p)
     conditionalFormatting(wb, sheet_name, cols = p_cols[2], rows = 2:(nrow(df) + 1), rule = "<0.05", style = style_p)
    
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
    #for (col in other_cols[other_cols %% 2 == 0]) {
    #    conditionalFormatting(wb, sheet_name, style = style_other, rule = "!=0", rows = 2:(nrow(df) + 1), cols = col, gridExpand = TRUE)
    #    conditionalFormatting(wb, sheet_name, style = style_other, rule = "==0", rows = 2:(nrow(df) + 1), cols = col, gridExpand = TRUE)
    #}

# Remove gridlines
showGridLines(wb, sheet_name, showGridLines = FALSE)
}
saveWorkbook(wb, paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_overlapping_geneFeatures/mto1_TEs_overlapping_geneFeatures.xlsx"), overwrite = T)
