library(dplyr)
library(openxlsx)

features_df <- data.frame(Feature = c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs"))



feature_count <- function(context=NULL, all = F) {
    x <- c()
    for (i.feature in c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs")) {
        if (!all) {
            x <- c(
                x,
                read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/", i.feature, "_", context, "_genom_annotations.csv")) %>%
                    # distinct(., gene_id) %>%  # count by DMRs. distinct for count by genes
                    nrow()
            )
        } else if (is.null(context)) {
            x <- c(
                x,
                rbind(
                    read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/", i.feature, "_CG_genom_annotations.csv")),
                    read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/", i.feature, "_CHG_genom_annotations.csv")),
                    read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/", i.feature, "_CHH_genom_annotations.csv"))
                ) %>%
                    # distinct(., gene_id) %>%  # count by DMRs. distinct for count by genes
                    nrow()
            )
        }
    }
    return(x)
}


features_df$CG = feature_count("CG")
features_df$CHG = feature_count("CHG")
features_df$CHH = feature_count("CHH")
features_df$All_contexts = feature_count(all=T)

features_df$Feature = gsub("fiveUTRs", "5'UTRs", features_df$Feature)
features_df$Feature = gsub("threeUTRs", "3'UTRs", features_df$Feature)

xl_headers <- names(features_df)
################
# save and edit EXCEL
wb <- createWorkbook()
# Define styles
cell_n_font_style <- createStyle(fontName = "Times New Roman", borderColour = "black")
last_row_style <- createStyle(fontName = "Times New Roman", border = "Bottom", borderStyle = "thick", borderColour = "black")
header_style <- createStyle(fontName = "Times New Roman", textDecoration = "bold", border = "Bottom", borderStyle = "thick", borderColour = "black")

sheet_name <- "Gene_features_count"
df <- features_df

addWorksheet(wb, sheet_name)
writeData(wb, sheet_name, df)

addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
addStyle(wb, sheet_name, style = last_row_style, rows = nrow(df) + 1, cols = 1:ncol(df), gridExpand = TRUE)
addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)

# Remove gridlines
showGridLines(wb, sheet_name, showGridLines = FALSE)

saveWorkbook(wb, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/Gene_features_count_by_DMRs.xlsx", overwrite = T)


