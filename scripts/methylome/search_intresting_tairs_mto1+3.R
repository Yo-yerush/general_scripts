#ann_vec = c("Genes","Promoters","CDS","Introns","fiveUTRs","threeUTRs","TEG")

#search_vec = c("methyltransferase","methylation","methionine")
#search_vec = "2.1.1."
#search_vec = tair_res_vec
#search_vec = "amino acid"
search_vec = "\\[GO:0016413\\]"
search_vec_merge = "O-acetyltransferase activity"

ann_vec = c("Genes","Promoters")

#ann_vec_for_plot = data.frame()
#lapply(year, function(x){ann_vec_for_plot[[x]][["year"]] <<- x})

if (length(grep("AT[0-9]G",search_vec[1])) == 0) {
  is.tair.list = F
  #search_vec_merge = ""
} else {
  is.tair.list = T
}

if (is.tair.list) {
  file_name_ending = paste0(paste(search_vec_merge, collapse = "_"),"_tair_list")
} else {
  file_name_ending = paste(search_vec, collapse = "_")
}

#tair_list = NA
tair_list_df_0 = data.frame(tair = NA, search_type = NA)
for (search_vec_loop in search_vec) {
  for (treat in c("mto1_vs_wt","mto3_vs_wt")) {
    for (cntx in c("CG","CHG","CHH")) {
      for (ann in ann_vec) {
        try({
          DMR_file = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_161123/",treat,"/genome_annotation/",cntx,"/",ann,"_",cntx,"_genom_annotations.csv"))
          
          DMR_intersting = DMR_file[with(DMR_file, grepl(search_vec_loop,#paste(search_vec, collapse = "|"), 
                                                         paste(
                                                           locus_tag, seqnames, start, end, width, strand, context, gene, pValue, regionType, db_xref, source, type, gbkey, gene_biotype, ID, note, gene_model_type, short_description, Curator_summary, Computational_description, EC.number, Gene.Ontology..biological.process., Gene.Ontology..molecular.function., Gene.Ontology..cellular.component., Tissue.specificity, Developmental.stage, Induction, Pathway, Function..CC., DOI.ID, PubMed.ID, sumReadsM1, sumReadsN1, proportion1, sumReadsM2, sumReadsN2, proportion2, cytosinesCount
                                                         ))),]
          
          #DMR_intersting = rbind(
          #  DMR_intersting,
          #  DMR_file[grep("2.1.1.",DMR_file$EC.number),]
          #)
          #tair_list = append(tair_list, unique(DMR_intersting$locus_tag))
          
          row_formula = (nrow(tair_list_df_0) + 1) : (nrow(tair_list_df_0) + length(unique(DMR_intersting$locus_tag)))
          tair_list_df_0[row_formula, 1] = unique(DMR_intersting$locus_tag)
          tair_list_df_0[row_formula, 2] = rep(search_vec_loop, length(unique(DMR_intersting$locus_tag)))
        })
      }
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
    for (treat in c("mto1","mto3")) {
      vec_of_names = append(vec_of_names, paste(treat,cntx,ann, sep = "_"))
    }
  }
}
vec_of_names = vec_of_names[-1]

tair_list_final = data.frame(matrix(ncol = length(vec_of_names)+2, nrow = nrow(tair_list_df)))
colnames(tair_list_final) = c("locus_tag","search_type",vec_of_names)
tair_list_final$locus_tag = tair_list_df$tair
tair_list_final$search_type = tair_list_df$search_type

j=3
while (j <= ncol(tair_list_final)) {
  split_col_name = stringr::str_split_1(names(tair_list_final)[j], "_")
  treat = split_col_name[1]
  cntx = split_col_name[2]
  ann = split_col_name[3]
  DMR_file = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_161123/",treat,"_vs_wt/genome_annotation/",cntx,"/",ann,"_",cntx,"_genom_annotations.csv"))
  
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
ar11 = read.csv("P:/yonatan/methionine/methylome_23/BSseq_pipeline/DMRcaller/Methylome.At/Methylome.At_scripts/data/Araport11_functional_descriptions_051023_short.txt", sep = "\t")
names(ar11)[1] = "locus_tag"
ar11 = ar11[!duplicated(ar11$locus_tag), grep("locus_tag|short_description",names(ar11))]
tair_and_description = merge.data.frame(ar11, heatmap_df, by = "locus_tag")
row.names(tair_and_description) = paste(tair_and_description$locus_tag, tair_and_description$short_description, sep = ": ")
tair_and_description = tair_and_description[order(tair_and_description$search_type),]
heatmap_mat = as.matrix(tair_and_description[,-c(1:3)])

if (length(grep("GO:",search_vec[1])) == 0) {
  file_name_ending = paste(search_vec, collapse = "_")
  labels_if_go = unique(tair_and_description$search_type)
  #search_vec_merge = ""
} else {
  file_name_ending = paste0(paste(search_vec_merge, collapse = "_"),"_GO_list")
  labels_if_go = paste(search_vec_merge, collapse = "_")
}

col_vec = c(rainbow(1, start = 0.1), rainbow(1, start = 0.43), rainbow(length(unique(tair_and_description$search_type)), start = 0.77, end = 0.92))
height = nrow(heatmap_mat)*0.25
if (height<4) {height=4}
library("ComplexHeatmap")
if (is.tair.list) {
  left_annotation = NULL
} else {
  left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = col_vec[3:(length(unique(tair_and_description$search_type))+2)]), 
                                 labels = labels_if_go))
}
pdf(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_161123/intresting_genes_",file_name_ending,".pdf"), width = 20, height = height)
Heatmap(heatmap_mat, na_col = "green", cluster_columns = F, cluster_column_slices = F, show_row_dend = F, show_heatmap_legend = F, 
        column_order = colnames(heatmap_mat), row_dend_reorder = F,
        row_names_max_width = unit(20, "cm"),
        rect_gp = gpar(col = "white", lwd = 2),
        column_gap = unit(3, "mm"), border = TRUE,#row_km = 3,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(rep(col_vec[1],3),rep(col_vec[2],3))),
                                                            labels = c("CG\nGenes", "CHG\nGenes", "CHH\nGenes",
                                                                       "CG\nPromoters", "CHG\nPromoters", "CHH\nPromoters"))),
        left_annotation = left_annotation,
        column_split = rep(1:6, each = 2),
        column_title=NULL,
        #(nrow(heatmap_mat)/3)
        row_split = tair_and_description$search_type,
        row_title=NULL,
        row_order = 1:nrow(heatmap_mat)
        #row_split = rep(1:7, each = 34)
        #column_km = 6
)
dev.off()
