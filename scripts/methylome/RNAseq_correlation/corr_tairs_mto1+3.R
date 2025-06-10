search_vec = "methyltransferase"
#search_vec = "histone"
#search_vec = "methylation"
#search_vec = "2.1.1.*"
#search_vec = "2.1.1.3[5-7][0-9]"
#search_vec = "DNA methyltransferase"
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
    DMR_G_or_L_0 = data.frame(treatment = NA, annotation = NA, context = NA, locus_tag = NA, regionType = NA)
    for (search_vec_loop in search_vec) {
      for (treat in c("mto1_vs_wt","mto3_vs_wt")) {
        for (cntx in c("CG","CHG","CHH")) {
          for (ann in c("Genes","Promoters")) {
            try({
              DMR_file = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treat,"/genome_annotation/",cntx,"/",ann,"_",cntx,"_genom_annotations.csv"))
              
              DMR_intersting = DMR_file[with(DMR_file, grepl(search_vec_loop,#paste(search_vec, collapse = "|"), 
                                                             paste(
                                                              EC.number, Gene.Ontology..biological.process., Gene.Ontology..molecular.function., Gene.Ontology..cellular.component.,regionType
                                                              #locus_tag, seqnames, start, end, width, strand, context, gene, pValue, regionType, db_xref, source, type, gbkey, gene_biotype, ID, note, gene_model_type, short_description, Curator_summary, Computational_description, EC.number, Gene.Ontology..biological.process., Gene.Ontology..molecular.function., Gene.Ontology..cellular.component., Tissue.specificity, Developmental.stage, Induction, Pathway, Function..CC., DOI.ID, PubMed.ID, sumReadsM1, sumReadsN1, proportion1, sumReadsM2, sumReadsN2, proportion2, cytosinesCount
                                                             ))),]
              row_formula = (nrow(tair_list_df_0) + 1) : (nrow(tair_list_df_0) + length(unique(DMR_intersting$locus_tag)))
              tair_list_df_0[row_formula, 1] = unique(DMR_intersting$locus_tag)
              tair_list_df_0[row_formula, 2] = rep(search_vec_loop, length(unique(DMR_intersting$locus_tag)))
              
              DMR_G_or_L_0 = rbind(DMR_G_or_L_0,data.frame(treatment = rep(treat, nrow(DMR_intersting)),
                                                           annotation = rep(ann, nrow(DMR_intersting)),
                                                           context = DMR_intersting$context,
                                                           locus_tag = DMR_intersting$locus_tag,
                                                           regionType = DMR_intersting$regionType))
            })
          }
        }
      }
    }
    tair_list = data.frame(locus_tag = unique(tair_list_df_0[-1,]$locus_tag))
    heatmap_df = merge.data.frame(tair_list,heatmap_df.0, by = "locus_tag")
    
    ncol_G_or_L = ncol(DMR_G_or_L_0)
    for (tmp.n in 1:nrow(DMR_G_or_L_0)) {
      DMR_G_or_L_0$tmp[tmp.n] = paste(DMR_G_or_L_0[tmp.n,1:ncol_G_or_L], collapse = "X")
    }
    DMR_G_or_L = DMR_G_or_L_0[!duplicated(DMR_G_or_L_0$tmp),]
    DMR_G_or_L = DMR_G_or_L[-1,-grep("tmp", names(DMR_G_or_L))]
  } else {
    heatmap_df = heatmap_df.0
  }
  
  
  ar11 = read.csv("P:/yonatan/methionine/methylome_23/BSseq_pipeline/DMRcaller/Methylome.At/Methylome.At_scripts/data/Araport11_functional_descriptions_051023_short.txt", sep = "\t")
  names(ar11)[1] = "locus_tag"
  ar11 = ar11[!duplicated(ar11$locus_tag), grep("locus_tag|short_description",names(ar11))]
  tair_and_description = merge.data.frame(ar11, heatmap_df, by = "locus_tag", all.y = T)
  

  
  tair_and_description[is.na(tair_and_description)] = 0
  tair_and_description_mto1 = tair_and_description[,grep("mto1",names(tair_and_description))]
  tair_and_description_mto3 = tair_and_description[,grep("mto3",names(tair_and_description))]
  
  tair_and_description_mto1 = tair_and_description_mto1[!apply(tair_and_description_mto1, 1, function(row) sum(row) == 0),]
  tair_and_description_mto3 = tair_and_description_mto3[!apply(tair_and_description_mto3, 1, function(row) sum(row) == 0),]
  
  
  row.names(tair_and_description) = paste(
    tair_and_description$locus_tag, tair_and_description$gene, tair_and_description$short_description, tair_and_description$transcript_id,
    sep = ": ")
  tair_and_description = tair_and_description[,-(1:4)]
  
  heatmap_mat_mto1 = as.matrix(tair_and_description_mto1)
  height_mto1 = nrow(heatmap_mat_mto1)*0.25
  if (height_mto1<4) {height_mto1=4}
  
  heatmap_mat_mto3 = as.matrix(tair_and_description_mto3)
  height_mto3 = nrow(heatmap_mat_mto3)*0.25
  if (height_mto3<4) {height_mto3=4}
  
  library("ComplexHeatmap")
  col_vec = c(rainbow(1, start = 0.1), rainbow(1, start = 0.43), rainbow(length(search_vec), start = 0.77, end = 0.92))
  if (!is.null(search_vec)) {
    left_annotation = rowAnnotation(foo = anno_block(#gp = gpar(fill = col_vec[3:(length(unique(tair_and_description$search_type))+2)]), 
                                                     labels = search_vec))
    search_vec_4file = gsub("\\*","all",search_vec)
    search_vec_4file = gsub("\\[.*\\]","num",search_vec_4file)
  } else {
    left_annotation = NULL
    search_vec_4file = ""
  }
  
  pdf(paste0("P:/yonatan/methionine/rnaseq_23/corr_with_methylations/mto_meth_RNA_corr_",search_vec_4file,".pdf"), width = 20, height = height_mto1+height_mto1+1)
  h1 = Heatmap(heatmap_mat_mto1, 
          #na_col = "green", 
          cluster_columns = FALSE, 
          cluster_rows = FALSE,
          cluster_column_slices = FALSE, 
          show_row_dend = FALSE, 
          show_heatmap_legend = FALSE, 
          column_order = colnames(heatmap_mat_mto1), 
          row_names_max_width = unit(20, "cm"),
          rect_gp = gpar(col = "white", lwd = 2),
          column_gap = unit(3, "mm"), 
          border = TRUE,
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(rep(col_vec[1],3),rep(col_vec[2],3))),
                                                              labels = c("CG\nGenes", "CHG\nGenes", "CHH\nGenes",
                                                                         "CG\nPromoters", "CHG\nPromoters", "CHH\nPromoters"))),
          left_annotation = left_annotation,
          column_split = 1:6,
          column_title = NULL,
          row_title = NULL,
          col = colorRamp2(c(min(heatmap_mat_mto3), 0, max(heatmap_mat_mto3)), c("blue", "gray95", "red"))
  )
  
  h2 = Heatmap(heatmap_mat_mto3, 
          #na_col = "green", 
          cluster_columns = FALSE, 
          cluster_rows = FALSE,
          cluster_column_slices = FALSE, 
          show_row_dend = FALSE, 
          show_heatmap_legend = FALSE, 
          column_order = colnames(heatmap_mat_mto3), 
          row_names_max_width = unit(20, "cm"),
          rect_gp = gpar(col = "white", lwd = 2),
          column_gap = unit(3, "mm"), 
          border = TRUE,
          #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(rep(col_vec[1],3),rep(col_vec[2],3))),
          #                                                    labels = c("CG\nGenes", "CHG\nGenes", "CHH\nGenes",
          #                                                               "CG\nPromoters", "CHG\nPromoters", "CHH\nPromoters"))),
          left_annotation = left_annotation,
          column_split = 1:6,
          column_title = NULL,
          row_title = NULL,
          col = colorRamp2(c(min(heatmap_mat_mto3), 0, max(heatmap_mat_mto3)), c("orange", "gray95", "green4"))
  )
  draw(h1 %v% h2)
  dev.off()
