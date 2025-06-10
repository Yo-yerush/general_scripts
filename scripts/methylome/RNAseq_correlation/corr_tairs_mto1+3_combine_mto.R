search_vec = "methyltransferase"
#search_vec = "methylation"
#search_vec = "methionine"

  vec_of_names = NA
  for (ann in c("genes","promoters")) {
    for (cntx in c("CG","CHG","CHH")) {
      for (treat in c("mto1","mto3")) {
        vec_of_names = append(vec_of_names, paste(treat,cntx,ann, sep = "_"))
      }
    }
  }
  vec_of_names = vec_of_names[-1]
  
  final_df = data.frame(transcript_id = NA)
  des_tmp = data.frame(transcript_id = NA, locus_tag = NA, gene = NA)
  for (name.l in vec_of_names) {
    split_col_name = stringr::str_split_1(name.l, "_")
    treat = split_col_name[1]
    cntx = split_col_name[2]
    ann = split_col_name[3]
    
    corr_file.0 = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/corr_with_methylations/",treat,"/",cntx,"/",ann,".corr.",cntx,".",treat,".csv"))
    
    des_tmp = rbind(des_tmp, corr_file.0[,grep("transcript_id|locus_tag|gene",names(corr_file.0))])
    #corr_file = corr_file.0[corr_file.0$pval < 0.05, grep("transcript_id|cor",names(corr_file.0))]
    corr_file.0 = corr_file.0[,grep("transcript_id|cor",names(corr_file.0))]
    corr_file_up= corr_file.0[corr_file.0$cor > 0.6,]
    corr_file_down = corr_file.0[corr_file.0$cor < -0.6,]
    corr_file = rbind(corr_file_up,corr_file_down)
    
    names(corr_file)[2] = name.l
    
    final_df = merge.data.frame(final_df,corr_file, by = "transcript_id", all = T)
  }
  des_tmp = des_tmp[-1,]
  des_tmp = des_tmp[!duplicated(des_tmp$transcript_id),]
  
  heatmap_df.0 = merge.data.frame(des_tmp,final_df, by = "transcript_id")
  
  if (!is.null(search_vec)) {
    tair_list_df_0 = data.frame(locus_tag = NA, search_type = NA)
    for (search_vec_loop in search_vec) {
      for (treat in c("mto1_vs_wt","mto3_vs_wt")) {
        for (cntx in c("CG","CHG","CHH")) {
          for (ann in c("Genes","Promoters")) {
            try({
              DMR_file = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treat,"/genome_annotation/",cntx,"/",ann,"_",cntx,"_genom_annotations.csv"))
              
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
    }
    tair_list = data.frame(locus_tag = unique(tair_list_df_0[-1,]$locus_tag))
    heatmap_df = merge.data.frame(tair_list,heatmap_df.0, by = "locus_tag")
  } else {
    heatmap_df = heatmap_df.0
  }
  
  
  ar11 = read.csv("P:/yonatan/methionine/methylome_23/BSseq_pipeline/DMRcaller/Methylome.At/Methylome.At_scripts/data/Araport11_functional_descriptions_051023_short.txt", sep = "\t")
  names(ar11)[1] = "locus_tag"
  ar11 = ar11[!duplicated(ar11$locus_tag), grep("locus_tag|short_description",names(ar11))]
  tair_and_description = merge.data.frame(ar11, heatmap_df, by = "locus_tag")
  
  row.names(tair_and_description) = paste(
    tair_and_description$locus_tag, tair_and_description$gene, tair_and_description$short_description, tair_and_description$transcript_id,
    sep = ": ")
  tair_and_description = tair_and_description[,-(1:4)]
  
  #tair_and_description = tair_and_description[order(tair_and_description$search_type),]
  
  
  heatmap_mat = as.matrix(tair_and_description)
  
  heatmap_mat[is.na(heatmap_mat)] <- 0
  
  col_vec = c(rainbow(1, start = 0.1), rainbow(1, start = 0.43), rainbow(length(search_vec), start = 0.77, end = 0.92))
  height = nrow(heatmap_mat)*0.25
  if (height<4) {height=4}
  
  library("ComplexHeatmap")
  
  if (!is.null(search_vec)) {
    left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = col_vec[3:(length(unique(tair_and_description$search_type))+2)]), 
                                                     labels = search_vec))
    search_vec_4file = search_vec
  } else {
    left_annotation = NULL
    search_vec_4file = ""
  }
  
  pdf(paste0("P:/yonatan/methionine/rnaseq_23/corr_with_methylations/mto_meth_RNA_corr_",search_vec_4file,".pdf"), width = 20, height = height)
  Heatmap(heatmap_mat, 
          #na_col = "green", 
          cluster_columns = FALSE, 
          cluster_rows = T,
          cluster_column_slices = FALSE, 
          show_row_dend = FALSE, 
          show_heatmap_legend = F, 
          column_order = colnames(heatmap_mat), 
          row_names_max_width = unit(20, "cm"),
          rect_gp = gpar(col = "white", lwd = 2),
          column_gap = unit(3, "mm"), 
          border = TRUE,
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(rep(col_vec[1],3),rep(col_vec[2],3))),
                                                              labels = c("CG\nGenes", "CHG\nGenes", "CHH\nGenes",
                                                                         "CG\nPromoters", "CHG\nPromoters", "CHH\nPromoters"))),
          left_annotation = left_annotation,
          column_split = rep(1:6, each = 2),
          column_title = NULL,
          row_title = NULL
  )
  
  dev.off()
