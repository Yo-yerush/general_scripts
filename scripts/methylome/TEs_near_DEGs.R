library(dplyr)
library(GenomicRanges)
library(openxlsx)

################################

log2FC_lim = 0.6

################################
###### TEs results
TE <- rbind(
  read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/Transposable_Elements_CG_genom_annotations.csv"),
  read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/Transposable_Elements_CHG_genom_annotations.csv"),
  read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/Transposable_Elements_CHH_genom_annotations.csv")
) %>%
  # dplyr::rename(gene_id = Transposon_Name) %>%
  mutate(tmp = paste(.$seqnames, .$start, .$end, .$gene_id, sep = "_del_")) %>%
  distinct(tmp, .keep_all = T) %>%
  dplyr::select(-tmp) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  sort()


################################
###### DEGs
DEGs = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv")
names(DEGs) = gsub("locus_tag", "gene_id", names(DEGs))
symbol_indx = DEGs %>% 
  select(gene_id, gene, log2FoldChange, pValue, short_description) %>%
  dplyr::rename(Symbol = gene)
symbol_indx$log2FoldChange = round(symbol_indx$log2FoldChange, 3)

upregulated = DEGs %>% 
  filter(log2FoldChange > log2FC_lim & padj < 0.05 & gene_model_type == "protein_coding") %>%
  select(gene_id)

downregulated = DEGs %>% 
  filter(log2FoldChange < -log2FC_lim & padj < 0.05 & gene_model_type == "protein_coding") %>%
  select(gene_id)

non_sig = DEGs %>% 
  filter(padj > 0.5 & gene_model_type == "protein_coding") %>%
  select(gene_id)


################################
###### TAIR10 annotations
gff3 = rtracklayer::import.gff3("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 gff3/TAIR10_GFF3_genes.gff")
gff3_df = gff3 %>%
  as.data.frame() %>%
  dplyr::rename(gene_id = ID) %>%
  filter(type == "gene") %>%
  select(seqnames,start,end,width,strand,gene_id)


################################
################################
###### overlapp with TEs up/down-stream

# overlap with TEs function
overlap_fun <- function(x.gr, x.TE=TE) {
  m = findOverlaps(x.gr, x.TE)
  x.out = x.gr[queryHits(m)] %>%
    as.data.frame() %>%
    mutate(tmp = paste(.$seqnames,.$start,.$end,.$strand,.$gene_id, sep = "_")) %>%
    distinct(tmp, .keep_all = T) %>%
    select(gene_id)
  
  # each gene and its TEs overlapped family
  TE2gene_df = data.frame(gene_id = x.gr[queryHits(m)]$gene_id,
                          TE_family = x.TE[subjectHits(m)]$Transposon_Family,
                          TE_super_family = x.TE[subjectHits(m)]$Transposon_Super_Family) %>%
    mutate(tmp = paste(.$gene_id, .$TE_family, .$TE_super_family, sep = "_del_")) %>%
    distinct(tmp, .keep_all = T) %>%
    select(-tmp)
  
  return(list(overlapIDs = x.out,
              TE2gene = TE2gene_df))
}

# remove suplicates fanilies or super-families in columns
remove_dup <- function(y) {
  y <- as.character(unique(unlist(strsplit(y, "; "))))
  paste(y, collapse = "; ")
}

#######
## main function
DEG_near_TEs <- function(deg_type, TE.upstream=F, TE.downstream=F, TE.both_sides=F) {
  deg.gr = merge(deg_type, gff3_df) %>% relocate(gene_id, .after = last_col()) %>% makeGRangesFromDataFrame(., keep.extra.columns = T)
  upstream = shift(deg.gr, 4000)
  downstream = shift(deg.gr, -4000)
  
  up_overlap = overlap_fun(upstream)$overlapIDs
  down_overlap = overlap_fun(downstream)$overlapIDs
  
  upstream_TE_fam =  overlap_fun(upstream)$TE2gene
  downstream_TE_fam =  overlap_fun(downstream)$TE2gene
  names(upstream_TE_fam) = gsub("TE", "upstream", names(upstream_TE_fam))
  names(downstream_TE_fam) = gsub("TE", "downstream", names(downstream_TE_fam))
  
  if (TE.both_sides) {
    return_df <- merge(up_overlap, down_overlap) %>%
      merge(., upstream_TE_fam) %>%
      merge(., downstream_TE_fam)
  } else if (TE.upstream) {
    return_df <- up_overlap %>%
      merge(., upstream_TE_fam)
  } else if (TE.downstream) {
    return_df <- down_overlap %>%
      merge(., downstream_TE_fam)
  }

  # group families and super-families by 'gene_id'
  return_df = return_df %>%
    group_by(gene_id) %>%
    summarise(
      across(contains("stream_family") | contains("stream_super_family"),
             ~remove_dup(paste(., collapse = "; ")))
    ) %>%
    merge(symbol_indx, ., by = "gene_id") %>%
    as.data.frame() %>%
    arrange(pValue)
  
  return(return_df)
}

xl_list <- list(
  Up.Stream = rbind(
    DEG_near_TEs(upregulated, TE.upstream = T),
    DEG_near_TEs(downregulated, TE.upstream = T)
  ),
  Down.Stream = rbind(
    DEG_near_TEs(upregulated, TE.downstream = T),
    DEG_near_TEs(downregulated, TE.downstream = T)
  ),
  Both.Sides = rbind(
    DEG_near_TEs(upregulated, TE.both_sides = T),
    DEG_near_TEs(downregulated, TE.both_sides = T)
  )
)

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

for (sheet_name in names(xl_list)) {
  addWorksheet(wb, sheet_name)
  
  xl_headers = names(xl_list[[sheet_name]])
  numeric_cols = grep("log2FoldChange|pValue", xl_headers)
  p_cols = grep("pValue", xl_headers)
  lfc_cols = grep("log2FoldChange", xl_headers)
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
  
}
saveWorkbook(wb, paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_near_DEGs/TEs_near_DEGs_201124.xlsx"), overwrite = T)
