library(dplyr)
library(openxlsx)

for (TE_near_DEG_type in c("Up.Stream","Both.Sides")) {
  
  TEs_upstream = read.xlsx("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_near_DEGs/TEs_near_DEGs_290924.xlsx",
                           sheet = TE_near_DEG_type)
  
  cor_bind <- function(cntx) {
    rbind(
      read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/mto1/",cntx,"/genes.corr.",cntx,".mto1.csv")) %>%
        mutate(cor.type = "genes") %>% select(-padj) %>% rename(gene_id = "locus_tag", cor.pval = "pval"),
      
      read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/corr_with_methylations/by_DEseq2/mto1/",cntx,"/promoters.corr.",cntx,".mto1.csv")) %>%
        mutate(cor.type = "promoters") %>% select(-padj) %>% rename(gene_id = "locus_tag", cor.pval = "pval")
    ) %>%
      merge.data.frame(., TEs_upstream, by = "gene_id") %>%
      mutate(Symbol = ifelse(is.na(Symbol), "", Symbol)) %>%
      arrange(pValue)
    
  }
  
  xl_list = list(CG = cor_bind("CG"),
                 CHG = cor_bind("CHG"),
                 CHH = cor_bind("CHH"))
  
  
  # save as excel
  clean_ASCII <- function(x) {
    x = gsub("\002", " ", x)
    x = gsub("\036", " ", x)
    
    #  x = gsub("[[:punct:]]", " ", x)
    #x = iconv(x, from = 'UTF-8', to = 'ASCII')
    return(x)
  }
  ################
  
  # save and edit EXCEL
  wb <- createWorkbook()
  modifyBaseFont(wb, fontName = "Times New Roman", fontSize = 11)
  
  for (sheet_name in names(xl_list)) {
    addWorksheet(wb, sheet_name, gridLines = F)
    
    xl_headers = names(xl_list[[sheet_name]])
    numeric_cols = grep("log2|pval|^cor$", tolower(xl_headers))
    p_cols = grep("pval", tolower(xl_headers))
    lfc_cols = grep("log2FoldChange|^cor$", xl_headers)
    fam_cols = grep("stream_family", xl_headers)
    super_fam_cols = grep("stream_super_family", xl_headers)
    ################
    
    # Define styles
    style_up <- createStyle(bgFill = "#f59d98")
    style_down <- createStyle(bgFill = "#c3ccf7")
    style_p <- createStyle(bgFill = "#f7deb0")
    style_fam <- createStyle(bgFill = "#e6ddf0")
    style_super_fam <- createStyle(bgFill = "#e1f5df")
    cell_n_font_style <- createStyle(border = "TopBottomLeftRight", borderColour = "black")
    header_style <- createStyle(textDecoration = "bold", border = "TopBottomLeftRight", borderStyle = "double")
    
    df <- xl_list[[sheet_name]]
    df <- data.frame(lapply(df, clean_ASCII))
    
    #### for make this columns as numeric
    df[,numeric_cols] = sapply(df[,numeric_cols], as.numeric)
    #### 
    
    writeData(wb, sheet_name, df)
    
    # headers and font
    addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
    addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
    
    # pValue
    conditionalFormatting(wb, sheet_name, cols = p_cols, rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
    
    # log2FC
    for (col in lfc_cols) {
      conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = ">0", style = style_up)
      conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = "<0", style = style_down)
    }
    
    # TEs groups style
    conditionalFormatting(wb, sheet_name, cols = fam_cols, rows = 2:(nrow(df)+1), rule = ">0", style = style_fam)
    conditionalFormatting(wb, sheet_name, cols = super_fam_cols, rows = 2:(nrow(df)+1), rule = ">0", style = style_super_fam)
    
    
    # Automatically set column widths based on the maximum content length in each column
    data.loop = xl_list[[sheet_name]]
    first.col = grep("short_description", names(data.loop))
    for (col in (first.col : ncol(data.loop))) {
      max_width <- max(c(nchar(as.character(data.loop[[col]])),
                         nchar(names(data.loop)[col])), na.rm = TRUE) + 1
      setColWidths(wb, sheet = sheet_name, cols = col, widths = max_width)
    }
  }
  saveWorkbook(wb, paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_near_DEGs/",TE_near_DEG_type,"_TE_n_cor_merged.xlsx"), overwrite = T)
}