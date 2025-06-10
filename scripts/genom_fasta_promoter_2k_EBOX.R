genom_fasta_by_refseq = function(id, # refseq id (XM or XP)
                                 id_as_LOC = F, # id as LOC and not refseq
                                 promoter = T,
                                 EBOX_motif = T) {
  
  if (id_as_LOC == T) {
    loc_id = id
  } else {
    loc_id = paste0("LOC",
                    rentrez::entrez_search("gene",id)$ids)
  }

  if (promoter == F) {
    EBOX_motif = F
  }
  
  date = gsub(" ","",format(Sys.time(), "%d %m %y"))
  time = gsub(" ","",format(Sys.time(), "%H %M"))
  date_time = paste(date, time, sep = "_")

  
  if (promoter == T) {
    file_name = "genom_seq_promoter_EBOX "
  } else if (promoter == F) {
    file_name = "genom_seq "
  }
  dir_name = paste0("P:/yonatan/saved.fasta.R/",id," ",file_name,date_time,"/")
  dir.create(dir_name)
  
  
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
if (promoter == T) {
  position = c(yo.1.gene_location[1]-2000,
               yo.1.gene_location[2]) + nchar(fasta_name_length)
} else {
  position = c(yo.1.gene_location[1],
                       yo.1.gene_location[2]) + nchar(fasta_name_length)
  }


yo.4 = substr(yo.3,position[1],position[2])

pattern = paste0(">",id," ",loc_id," ", yo.1.gene_name," ",fasta_ch_name,"\n")
yo.5 = paste0(pattern,yo.4)

yo.6 = substr(yo.5,
                   1+nchar(pattern),
                   2000+nchar(pattern))
motif_txt = c(rep("",4))
if (EBOX_motif == T) {
  motif_txt = c(
  try(if (grep("CAGCTG|CACCTG",yo.6) == 1) {
    motif_txt[1] = "A"
  }, silent = T),
  try(if (grep("CACGTG|CATGTTG",yo.6) == 1) {
    motif_txt[2] = "B"
  }, silent = T),
  try(if (grep("ACGTG|GCGTG",yo.6) == 1) {
    motif_txt[3] = "C"
  }, silent = T),
  try(if (grep("CACGCG|CACGAG",yo.6) == 1) {
    motif_txt[4] = "E"
  }, silent = T)
   )
}
motif_txt = gsub("Error(.*)","",motif_txt)
if (EBOX_motif == T) {
  motif_txt = gsub("Error(.*)","",motif_txt)
  if (nchar(motif_txt[1])+nchar(motif_txt[2])+nchar(motif_txt[3])+nchar(motif_txt[4]) != 0) {
    motif_txt = gsub("Error(.*)","",motif_txt)
    write(c("group",motif_txt),
                paste0(dir_name,id,"_",loc_id,"_EBOX_motif_class",".txt"))
  }

}
write(yo.5, paste0(dir_name,id,"_",loc_id,"_",file_name, ".fasta"))
}
