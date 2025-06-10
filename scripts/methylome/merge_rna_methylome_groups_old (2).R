library(dplyr)
library(writexl)
library(openxlsx)

###############
## fix it like mto1 paper
##############

xl_groups <- function(treatment) {
  
  all_res = as.data.frame(readxl::read_xlsx("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_genes.xlsx", sheet = treatment))
  
  xl_list = list(
    EC_211 = all_res[grep("^2\\.1\\.1\\.", tolower(all_res$EC.number)),] %>% filter(RNA_padj < 0.05),
    EC_21137 = all_res[grep("^2\\.1\\.1\\.37", tolower(all_res$EC.number)),],
    DNA_C5_MT = all_res[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(all_res$Protein.names)),],
    SAM_MT = all_res[grep("sam-binding methyltransferase", tolower(all_res$Protein.families)),] %>% filter(RNA_padj < 0.05),
    DNA_meth_function = all_res[grep("dna methylation", tolower(all_res$Function..CC.)),] %>% filter(RNA_padj < 0.05),
    DNA_meth_bp = all_res[grep("dna methylation", tolower(all_res$Gene.Ontology..biological.process.)),] %>% filter(RNA_padj < 0.05),
    chromosome_CC = all_res[grep("chromosome", tolower(all_res$Gene.Ontology..cellular.component.)),] %>% filter(RNA_padj < 0.05),
    chromatin_CC = all_res[grep("chromatin", tolower(all_res$Gene.Ontology..cellular.component.)),] %>% filter(RNA_padj < 0.05),
    nucleosome_CC = all_res[grep("nucleosome", tolower(all_res$Gene.Ontology..cellular.component.)),] %>% filter(RNA_padj < 0.05),
    chromatin_remodeling_BP = all_res[grep("chromatin remodeling", tolower(all_res$Gene.Ontology..biological.process.)),] %>% filter(RNA_padj < 0.05),
    chromatin_remodeling_pro_name = all_res[grep("chromatin remodeling", tolower(all_res$Protein.names)),] %>% filter(RNA_padj < 0.05),
    histone_H3_acetylation_BP = all_res[grep("histone H3 acetylation",tolower(all_res$Gene.Ontology..biological.process.)),] %>% filter(RNA_padj < 0.05),
    cell_wall_BP = all_res[grep("cell wall",tolower(all_res$Gene.Ontology..biological.process.)),] %>% filter(RNA_padj < 0.05),
    AA_biosynthesis = all_res[grep("amino-acid biosynthesis",tolower(all_res$Pathway)),] %>% filter(RNA_padj < 0.05),
    #AA_and_meth_biosynthesis = all_res[grep("amino-acid biosynthesis&methionine",tolower(all_res$Pathway)),] %>% filter(RNA_padj < 0.05),
    meth_biosynthesis = all_res[grep("methionine biosynthesis",tolower(all_res$Pathway)),] %>% filter(RNA_padj < 0.05),
    meth_bp = all_res[grep("methionine",tolower(all_res$Gene.Ontology..biological.process.)),] %>% filter(RNA_padj < 0.05),
    meth_mf = all_res[grep("methionine",tolower(all_res$Gene.Ontology..molecular.function.)),] %>% filter(RNA_padj < 0.05),
    aspartate_kinase_activity_MF = all_res[grep("aspartate kinase activity",tolower(all_res$Gene.Ontology..molecular.function.)),] %>% filter(RNA_padj < 0.05),
    cysteine_bp = all_res[grep("cysteine",tolower(all_res$Gene.Ontology..biological.process.)),] %>% filter(RNA_padj < 0.05),
    serine_bp = all_res[grep("serine",tolower(all_res$Gene.Ontology..biological.process.)),] %>% filter(RNA_padj < 0.05),
    GST_genes = all_res[grep("^GST",all_res$gene),],
    GSTU_genes = all_res[grep("^GSTU",all_res$gene),],
    sulfur_mf = all_res[grep("sulfur",tolower(all_res$Gene.Ontology..molecular.function.)),] %>% filter(RNA_padj < 0.05),
    sulfur_function = all_res[grep("sulfur",tolower(all_res$Function..CC.)),] %>% filter(RNA_padj < 0.05)
  )
  
  ########################################################
  
  ################
  clean_ASCII <- function(x) {
    x = gsub("\002", " ", x)
    x = gsub("\036", " ", x)
    
    #  x = gsub("[[:punct:]]", " ", x)
    #x = iconv(x, from = 'UTF-8', to = 'ASCII')
    return(x)
  }
  ################
  
  ################
  xl_headers = names(xl_list[[1]])
  numeric_cols = grep("RNA_|CG_|CHG_|CHH_", xl_headers)
  p_cols = grep("RNA_p", xl_headers)
  lfc_cols = grep("RNA_log2FC|CG_|CHG_|CHH_", xl_headers)
  other_cols = grep("Curator_summary", xl_headers) : length(xl_headers)
  ################
  
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
    #  setColWidths(wb, sheet_name, cols = other_cols, widths = 10)
    #  setColWidths(wb, sheet_name, cols = lfc_cols, widths = 4)
    #  setColWidths(wb, sheet_name, cols = p_cols, widths = 6)
    
    df <- xl_list[[sheet_name]]
    df <- data.frame(lapply(df, clean_ASCII))
    df[,numeric_cols] = sapply(df[, numeric_cols], as.numeric)
    
    writeData(wb, sheet_name, df)
    
    addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
    addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
    
    conditionalFormatting(wb, sheet_name, cols = p_cols[1], rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
    conditionalFormatting(wb, sheet_name, cols = p_cols[2], rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
    
    for (col in lfc_cols) {
      conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = ">0", style = style_up)
      conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = "<0", style = style_down)
    }
    
    for (col in other_cols[other_cols %% 2 == 1]) {
      conditionalFormatting(wb, sheet_name, style = style_other, rule = "!=0", rows = 2:(nrow(df)+1), cols = col, gridExpand = TRUE)
      conditionalFormatting(wb, sheet_name, style = style_other, rule = "==0", rows = 2:(nrow(df)+1), cols = col, gridExpand = TRUE)
    }
  }
  saveWorkbook(wb, paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_", treatment, "_groups.xlsx"), overwrite = T)
}

xl_groups("mto1")
xl_groups("mto3")
xl_groups("35s")
