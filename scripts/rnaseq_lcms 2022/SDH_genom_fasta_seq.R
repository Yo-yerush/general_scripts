genom_fasta_by_refseq = function(id, id_as_LOC = F, EBOX_motif = T) {
  
  if (id_as_LOC == T) {
    loc_id = id
  } else {
    loc_id = paste0("LOC",
                    rentrez::entrez_search("gene",id)$ids)
  }

  
  
yo.1 = rentrez::entrez_fetch(db = "gene", id = loc_id, 
                             rettype = "fasta")
yo.1.1 = gsub("[\r\n]","",yo.1)
yo.1.gene_name = stringr::str_match(yo.1.1, "Designations: \\s*(.*?)\\s*Chromosome:")[,2]
yo.1.chr_name = stringr::str_extract(yo.1.1, "NC_(\\d+)\\.1")
yo.1.gene_location = stringr::str_extract(yo.1.1, "(\\d+)\\.\\.(\\d+)")
yo.1.gene_location = strsplit(yo.1.gene_location, "\\.\\.")
yo.1.gene_location = as.integer(yo.1.gene_location[[1]])

yo.2 = rentrez::entrez_fetch(db = "nucleotide", id = yo.1.chr_name, 
                             rettype = "fasta")
yo.3 = gsub("[\r\n]","",yo.2)
fasta_name_length = stringr::str_extract(yo.3, ">NC_(\\d+)\\.1(.*)sequence")
fasta_ch_name = substr(yo.3,2,12)
position = c(yo.1.gene_location[1]-3000,
             yo.1.gene_location[2]) + nchar(fasta_name_length)

yo.4 = substr(yo.3,position[1],position[2])

pattern = paste0(">",id," ",loc_id," ", yo.1.gene_name," ",fasta_ch_name,"\n")
yo.5 = paste0(pattern,yo.4)

yo.6 = substr(yo.5,
                   1+nchar(pattern),
                   3000+nchar(pattern))
motif_txt = c(rep("",4))
if (EBOX_motif == T) {
  motif_txt = c(
  try(if (grep("CAGCTG|CACCTG",yo.6) == 1) {
    motif_txt[1] = "group A"
  }, silent = T),
  try(if (grep("CACGTG|CATGTTG",yo.6) == 1) {
    motif_txt[2] = "group B"
  }, silent = T),
  try(if (grep("ACGTG|GCGTG",yo.6) == 1) {
    motif_txt[3] = "group C"
  }, silent = T),
  try(if (grep("CACGCG|CACGAG",yo.6) == 1) {
    motif_txt[4] = "group E"
  }, silent = T)
   )
}
motif_txt = gsub("Error(.*)","",motif_txt)
if (EBOX_motif == T) {
  motif_txt = gsub("Error(.*)","",motif_txt)
  if (nchar(motif_txt[1])+nchar(motif_txt[2])+nchar(motif_txt[3])+nchar(motif_txt[4]) != 0) {
    motif_txt = gsub("Error(.*)","",motif_txt)
    write.table(motif_txt,paste0("P:/yonatan/RNAseq_yonatan_2021/SDH_variants_RNAseq_yonatan/SDH_genom_fasta/",id," ",loc_id," EBOX_motifs.txt"),
                quote = F, row.names = F)
  }

}
write.table(yo.5,paste0("P:/yonatan/RNAseq_yonatan_2021/SDH_variants_RNAseq_yonatan/SDH_genom_fasta/",id," ",loc_id,".fasta"),
            quote = F, row.names = F)
}
