library(dplyr)
library(writexl)
library(openxlsx)


  xl_groups <- function(treatment) {

  all_res = as.data.frame(readxl::read_xlsx("P:/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_transcripts_250123.xlsx", sheet = treatment))
  
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
  # save and edit EXCEL
  wb <- createWorkbook()
  # Define styles
  style_up <- createStyle(bgFill = "#f59d98", fontSize = 11, fontName = "Times New Roman")
  style_down <- createStyle(bgFill = "#c3ccf7", fontSize = 11, fontName = "Times New Roman")
  style_p <- createStyle(bgFill = "#f7deb0", fontSize = 11, fontName = "Times New Roman")
  style_other <- createStyle(bgFill = "#daf7d7", fontSize = 11, fontName = "Times New Roman")
  header_style <- createStyle(fontSize = 12, fontName = "Times New Roman", textDecoration = "bold", border = "Bottom", borderStyle = "double")
  cell_n_font_style <- createStyle(fontSize = 11, fontName = "Times New Roman", border = "TopBottomLeftRight", borderColour = "black")
  
  for (sheet_name in names(xl_list)) {
    addWorksheet(wb, sheet_name)
    df <- xl_list[[sheet_name]]
    writeData(wb, sheet_name, df)
    
    addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
    addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
    
    conditionalFormatting(wb, sheet_name, cols = 4, rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
    for (col in c(3, 5:10)) {
      conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = ">0", style = style_up)
      conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = "<0", style = style_down)
    }
    
    other_cols = 11:ncol(df)
    for (col in other_cols[other_cols %% 2 == 1]) {
      conditionalFormatting(wb, sheet_name, style = style_other, rule = "!=0", rows = 2:(nrow(df)+1), cols = col, gridExpand = TRUE)
      conditionalFormatting(wb, sheet_name, style = style_other, rule = "==0", rows = 2:(nrow(df)+1), cols = col, gridExpand = TRUE)
    }
  }
  saveWorkbook(wb, paste0("P:/yonatan/methionine/NGS_merged_results/merged_results_", treatment, "_groups.xlsx"), overwrite = TRUE)
  
  rm(xl_list)
}
  
  xl_groups("mto1")
  xl_groups("mto3")