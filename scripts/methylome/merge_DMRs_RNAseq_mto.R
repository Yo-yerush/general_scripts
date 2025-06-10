# first run the function script in the bottom of this file
search_vec_merge = "amino acid"

DMRs_mto1 = search_interesting_genes(search_vec_merge,"mto1_vs_wt")
DMRs_mto3 = search_interesting_genes(search_vec_merge,"mto3_vs_wt")

mto1_rnaseq = read.csv("P:/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv")
mto1_rnaseq = mto1_rnaseq[mto1_rnaseq$padj < 0.05,]

mto3_rnaseq = read.csv("P:/yonatan/methionine/rnaseq_23/met23/mto3_vs_wt/all.transcripts.mto3_vs_wt.DE.csv")
mto3_rnaseq = mto3_rnaseq[mto3_rnaseq$padj < 0.05,]

merged_mto1 = merge.data.frame(mto1_rnaseq, DMRs_mto1, by = "locus_tag")
merged_mto3 = merge.data.frame(mto3_rnaseq, DMRs_mto3, by = "locus_tag")

tair_res_vec = unique(c(merged_mto1$locus_tag, merged_mto3$locus_tag))
# can use 'tair_res_vec' to run in 'search_intresting_tairs_mto1+3'search_intresting_tairs_mto1+3 script
# as: 'search_vec = tair_res_vec'

merged_mto = merge.data.frame(merged_mto1, merged_mto3, by = "locus_tag")





search_interesting_genes <- function(search_vec,treat,
                                     Methylome.At_results = "P:/yonatan/methionine/methylome_23/BSseq_results_161123") {
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
  heatmap_df = heatmap_df[order(heatmap_df$search_type),]
 
  return(heatmap_df)
}


