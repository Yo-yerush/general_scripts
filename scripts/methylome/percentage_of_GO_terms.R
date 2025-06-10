#####
{
  library(dplyr)
  library(writexl)
  library(openxlsx)
  
  treatment = "mto1"
  all_res = as.data.frame(readxl::read_xlsx("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_genes.xlsx", sheet = treatment))
 # all_res = all_res[,-ncol(all_res)] %>%
#    relocate("transcript_id", .after = "locus_tag")
  
  ################
  offspring_fun <- function(go_id, xx = as.list(GO.db::GOBPOFFSPRING)) { # 'GOBPCHILDREN' for child terms
    
    child_terms_0 = as.character(xx[[go_id]])
    child_terms = child_terms_0
    
    for (i in 1:length(child_terms_0)) {
      child_terms = c(child_terms, as.character(xx[[child_terms[i]]]))
    } 
    
    return(child_terms[!is.na(child_terms)] %>% unique()) # %>% paste(collapse = "|"))
  }
  
  child_terms_response = offspring_fun("GO:0006950") # response to stress
  child_terms_stimulus = offspring_fun("GO:0050896") # response to stimulus
  child_terms_biotic = offspring_fun("GO:0009607") # response to biotic stimulus
  child_terms_abiotic = offspring_fun("GO:0009628") # response to abiotic stimulus
  child_terms_bio_p = offspring_fun("GO:0009058") # biosynthetic process
  child_terms_aa = offspring_fun("GO:0008652") # amino acid biosynthetic process
  child_terms_methionine = offspring_fun("GO:0009086") # methionine biosynthetic process
  child_terms_chromatin_org = offspring_fun("GO:0006325") # chromatin organization
  child_terms_chromatin_rem = offspring_fun("GO:0006338") # chromatin remodeling
  ################
  
  grep_position <- function(x) {
    vec = NULL
    for (terms_l in x) {
      vec = c(vec, grep(terms_l, all_res$Gene.Ontology..biological.process.))
    }
    return(unique(vec))
  }
  
  ################
  xl_list = list(
    response_to_stress = all_res[grep_position(child_terms_response),], # %>% filter(RNA_pvalue < 0.05) # 'response to stress' child terms
    response_to_stimulus = all_res[grep_position(child_terms_stimulus),],
    response_to_biotic_stimulus = all_res[grep_position(child_terms_biotic),],
    response_to_abiotic_stimulus = all_res[grep_position(child_terms_abiotic),],
    biosynthetic_process = all_res[grep_position(child_terms_bio_p),],
    amino_acid_biosynthetic_process = all_res[grep_position(child_terms_aa),],
    methionine_biosynthetic_process = all_res[grep_position(child_terms_methionine),],
    #chromatin_organization = all_res[grep_position(child_terms_chromatin_org),]
    chromatin_remodeling = all_res[grep_position(child_terms_chromatin_rem),] # %>% filter(RNA_pvalue < 0.05) # 'defense response' related terms
  )
}
#####


for (term_name in names(xl_list)) {
  a = xl_list[[term_name]] %>% filter(RNA_padj < 0.05) %>% nrow()
  b = all_res %>% nrow()
  c = all_res %>% filter(RNA_padj < 0.05) %>% nrow()
  (a/b)*100
  (a/c)*100
  
  # a - 'response_to_stress' sig. transcripts
  # b - all transcripts
  # c - all sig. transcripts
  
  
  
  d = xl_list[[term_name]] %>% filter(RNA_padj < 0.05) %>% distinct(locus_tag, .keep_all = T) %>% nrow()
  e = all_res %>% distinct(locus_tag, .keep_all = T) %>% nrow()
  f = all_res %>% distinct(locus_tag, .keep_all = T) %>% filter(RNA_padj < 0.05) %>% nrow()
  (d/e)*100
  (d/f)*100
  
  # d - 'response_to_stress' sig. genes
  # e - all genes
  # f - all sig. genes
  
  transcripts = data.frame(total = b,
                           sig = c,
                           term = a,
                           term.VS.total = paste0(round((a/b)*100,2),"%"),
                           term.VS.sig = paste0(round((a/c)*100,2),"%"))
  
  genes = data.frame(total = e,
                     sig = f,
                     term = d,
                     term.VS.total = paste0(round((d/e)*100,2),"%"),
                     term.VS.sig = paste0(round((d/f)*100,2),"%"))
  
  return_df = rbind(transcripts, genes)
  row.names(return_df) = c("transcripts", "genes")
  
  message(gsub("_"," ",term_name))

  #cat(term_name, "\n")
  print(return_df)
  cat("\n\n")
}
  
