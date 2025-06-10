search_interesting_genes <- function(search_vec,
                                     treat,
                                     Methylome.At_results = "./results",
                                     path_2_results = "./results",
                                     Araport11_description_file = "./Methylome.At_scripts/data/Araport11_functional_descriptions_051023_short.txt") {
  add_d <- function(x) {
    if (substr(x, nchar(x), nchar(x)) == "/") {
      x = substr(x, 1, nchar(x)-1)
    }
    return(x)
  }
  
  ann_vec = c("Genes","Promoters")
  
  tair_list_df_0 = data.frame(tair = NA, search_type = NA)
  for (search_vec_loop in search_vec) {
    for (cntx in c("CG","CHG","CHH")) {
      for (ann in ann_vec) {
        try({
          DMR_file = read.csv(paste0(Methylome.At_results,"/",treat,"/genome_annotation/",cntx,"/",ann,"_",cntx,"_genom_annotations.csv"))
          
          DMR_intersting = DMR_file[with(DMR_file, grepl(search_vec_loop,#paste(search_vec, collapse = "|"), 
                                                         paste(
                                                           locus_tag, seqnames, start, end, width, strand, context, gene, pValue, regionType, db_xref, source, type, gbkey, gene_biotype, ID, note, gene_model_type, short_description, Curator_summary, Computational_description, EC.number, Gene.Ontology..biological.process., Gene.Ontology..molecular.function., Gene.Ontology..cellular.component., Tissue.specificity, Developmental.stage, Induction, Pathway, Function..CC., DOI.ID, PubMed.ID, sumReadsM1, sumReadsN1, proportion1, sumReadsM2, sumReadsN2, proportion2, cytosinesCount
                                                         ))),]

          
          row_formula = (nrow(tair_list_df_0) + 1) : (nrow(tair_list_df_0) + length(unique(DMR_intersting$locus_tag)))
          tair_list_df_0[row_formula, 1] = unique(DMR_intersting$locus_tag)
          tair_list_df_0[row_formula, 2] = rep(search_vec_loop, length(unique(DMR_intersting$locus_tag)))
        })
      }
    }
  }
  tair_list_df_0 = tair_list_df_0[-1,]
  #tair_list_df_0 = tair_list_df_0[!duplicated(tair_list_df_0$tair),]
  
  if (length(search_vec) > 1) {
    tair_list_df = data.frame(tair = unique(tair_list_df_0$tair), search_type = NA)
    
    i=1
    while (i <= length(tair_list_df$tair)) {
      tair_list_df$search_type[i] = paste(unique(tair_list_df_0$search_type[grep(tair_list_df$tair[i], tair_list_df_0$tair)]),
                                          collapse = "+")
      i=i+1
    }
  } else {tair_list_df = tair_list_df_0[!duplicated(tair_list_df_0$tair),]}
  
  
  vec_of_names = NA
  for (ann in ann_vec) {
    for (cntx in c("CG","CHG","CHH")) {
      vec_of_names = append(vec_of_names, paste(treat,cntx,ann, sep = " "))
    }
  }
  vec_of_names = vec_of_names[-1]
  
  tair_list_final = data.frame(matrix(ncol = length(vec_of_names)+2, nrow = nrow(tair_list_df)))
  colnames(tair_list_final) = c("locus_tag","search_type",vec_of_names)
  tair_list_final$locus_tag = tair_list_df$tair
  tair_list_final$search_type = tair_list_df$search_type
  
  j=3
  while (j <= ncol(tair_list_final)) {
    split_col_name = stringr::str_split_1(names(tair_list_final)[j], " ")
    #treat = split_col_name[1]
    cntx = split_col_name[2]
    ann = split_col_name[3]
    DMR_file = read.csv(paste0(Methylome.At_results,"/",treat,"/genome_annotation/",cntx,"/",ann,"_",cntx,"_genom_annotations.csv"))
    
    i=1
    while (i <= nrow(tair_list_final)) {
      q1_g = grep(tair_list_final$locus_tag[i], DMR_file$locus_tag)
      q1 = unique(DMR_file$regionType[q1_g])
      if (length(q1) == 0) {
        tair_list_final[i,names(tair_list_final)[j]] = 0
      } else if (length(q1) == 1) {
        q2 = ifelse(q1 == "gain", 1, -1)
        tair_list_final[i,names(tair_list_final)[j]] = q2
      } else if (length(q1) == 2) {
        tair_list_final[i,names(tair_list_final)[j]] = NA
      }
      i=i+1 
    }
    j=j+1
  }
  
  heatmap_df = tair_list_final
  ar11 = read.csv(Araport11_description_file, sep = "\t")
  names(ar11)[1] = "locus_tag"
  ar11 = ar11[!duplicated(ar11$locus_tag), grep("locus_tag|short_description",names(ar11))]
  tair_and_description = merge.data.frame(ar11, heatmap_df, by = "locus_tag")
  row.names(tair_and_description) = paste(tair_and_description$locus_tag, tair_and_description$short_description, sep = ": ")
  tair_and_description = tair_and_description[order(tair_and_description$search_type),]
  heatmap_mat = as.matrix(tair_and_description[,-c(1:3)])
  
  col_vec = c(rainbow(1, start = 0.1), rainbow(1, start = 0.43), rainbow(length(unique(tair_and_description$search_type)), start = 0.77, end = 0.92))
  #col_vec = rainbow(length(unique(tair_and_description$search_type))+6, start = 0.42, end)
  #col_vec = col_vec[-c(1, 2, length(col_vec)-1, length(col_vec)-2)]
  library("ComplexHeatmap")
  
# change file name if is allready existing
  for (i in 1:1000) {
    if (file.exists(paste0(path_2_results,"/intresting_genes_",i,"_",treat,".pdf")) == F) {
      new_file_name = paste0(path_2_results,"/intresting_genes_",i,"_",treat,".pdf")
      break
    }
  }

  pdf(new_file_name, width = 14, height = (nrow(heatmap_mat)*0.25), family = "serif")
  Heatmap(heatmap_mat, na_col = "green", cluster_columns = F, cluster_column_slices = F, show_row_dend = F, show_heatmap_legend = F, 
          column_order = colnames(heatmap_mat), row_dend_reorder = F,
          row_names_max_width = unit(22, "cm"),
          rect_gp = gpar(col = "white", lwd = 2),
          column_gap = unit(3, "mm"), row_gap = unit(2, "mm"), 
          border = TRUE,#row_km = 3,
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(rep(col_vec[1],3),rep(col_vec[2],3))),
                                                              labels = c("CG\nGenes", "CHG\nGenes", "CHH\nGenes",
                                                                         "CG\nPromoters", "CHG\nPromoters", "CHH\nPromoters"))),
          left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = col_vec[3:(length(unique(tair_and_description$search_type))+2)]), 
                                                           labels = unique(tair_and_description$search_type))),
          column_split = 1:6,
          column_title = gt_render(paste0("<span style='font-size:12pt'>**",gsub("_"," ",treat),"**</span>"), 
                                   r = unit(2, "pt"), 
                                   padding = unit(c(2, 2, 2, 2), "pt")), 
          #column_title = gsub("_","",treat),
          
          #(nrow(heatmap_mat)/3)
          row_split = tair_and_description$search_type,
          row_title = NULL,
          row_order = 1:nrow(heatmap_mat)
  )
  dev.off()
  
}

search_interesting_genes(search_vec = c("methyltransferase","methylation","methionine"),
                         treat = "mto1_vs_wt",
                         Araport11_description_file = "P:/yonatan/methionine/methylome_23/BSseq_pipeline/DMRcaller/Methylome.At/Methylome.At_scripts/data/Araport11_functional_descriptions_051023_short.txt",
                         Methylome.At_results = "P:/yonatan/methionine/methylome_23/BSseq_results_161123",
                         path_2_results = "P:/yonatan/methionine/methylome_23/BSseq_results_161123")

search_interesting_genes(search_vec = c("methyltransferase","methylation","methionine"),
                         treat = "mto3_vs_wt",
                         Araport11_description_file = "P:/yonatan/methionine/methylome_23/BSseq_pipeline/DMRcaller/Methylome.At/Methylome.At_scripts/data/Araport11_functional_descriptions_051023_short.txt",
                         Methylome.At_results = "P:/yonatan/methionine/methylome_23/BSseq_results_161123",
                         path_2_results = "P:/yonatan/methionine/methylome_23/BSseq_results_161123")


search_vec = c("methyltransferase","methylation","methionine")
treat = "mto1_vs_wt"
Araport11_description_file = "P:/yonatan/methionine/methylome_23/BSseq_pipeline/DMRcaller/Methylome.At/Methylome.At_scripts/data/Araport11_functional_descriptions_051023_short.txt"
Methylome.At_results = "P:/yonatan/methionine/methylome_23/BSseq_results_161123"
path_2_results = "P:/yonatan/methionine/methylome_23/BSseq_results_161123"