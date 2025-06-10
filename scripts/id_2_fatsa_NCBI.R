id_2_fasta <- function(id,
                       input_seq = "mRNA or protein",
                       output_seq = "mRNA or protein",
                       save_FASTA = F,
                       full_name_file = T
) {

  library(rentrez)
  library(stringr)
  
  ##################################
  if (input_seq == "mRNA or protein") {
    stop("Must specify either (not both) 'mRNA' or 'protein' as 'input_seq'", call. = F)
  }
  if (output_seq == "mRNA or protein") {
    stop("Must specify either (not both) 'mRNA' or 'protein' as 'output_seq'", call. = F)
  }
  
  ##################################
  
  if (input_seq == "mRNA") {
    if (output_seq == "mRNA") {
      
      fasta_seq = entrez_fetch(db = "nuccore", id = id, 
                               rettype = "fasta")
      
    } else if (output_seq == "protein") {
      
      new_id = entrez_link(dbfrom = "nuccore",
                            id = id,
                            db = "protein")
      try({
        fasta_seq = entrez_fetch(db = "protein", id = new_id$links$nuccore_protein, 
                                   rettype = "fasta")},
        silent = T)
      try({
        fasta_seq = entrez_fetch(db = "protein", id = new_id[[1]]$links$nuccore_protein, 
                                                rettype = "fasta")},
        silent = T)
      
    } 
  }
  
  ##################################
  
  if (input_seq == "protein") {
    if (output_seq == "mRNA") {
      
      new_id = entrez_link(dbfrom = "protein",
                           id = id,
                           db = "nuccore")
      try({
        fasta_seq = entrez_fetch(db = "nuccore", id = new_id$links$protein_nuccore_mrna, 
                                   rettype = "fasta")},
        silent = T)
      try({
        fasta_seq = entrez_fetch(db = "nuccore", id = new_id[[1]]$links$protein_nuccore_mrna, 
                                                rettype = "fasta")},
        silent = T)
      
    } else if (output_seq == "protein") {
      
      fasta_seq = entrez_fetch(db = "protein", id = id, 
                               rettype = "fasta")
      
    } 
  }
  
  ##################################
  
  if (output_seq == "mRNA") {
    fasta_name = substr(fasta_seq, 2,nchar(str_extract(fasta_seq, ">(.*), mRNA")))
    if (full_name_file == T) {
      fasta_name = gsub(":|-|\\/","_",fasta_name)
    } else {fasta_name = substr(fasta_name,1,14)}
  } else if (output_seq == "protein") {
    fasta_name = substr(fasta_seq, 2,nchar(str_extract(fasta_seq, ">(.*)\\]")))
    if (full_name_file == T) {
      fasta_name = gsub(":|-|\\/", "_", fasta_name)
    } else {fasta_name = substr(fasta_name,1,14)}
  }
  
if (save_FASTA == F) {
  if (output_seq == "mRNA") {
    print_fasta = gsub("[\r\n]","",fasta_seq)
    print_fasta_name = substr(print_fasta, 1,nchar(str_extract(print_fasta, ">(.*), mRNA")))
    print_fasta_seq = substr(print_fasta, 1+nchar(str_extract(print_fasta, ">(.*), mRNA")),nchar(print_fasta))
  } else if (output_seq == "protein") {
    print_fasta = gsub("[\r\n]","",fasta_seq)
    print_fasta_name = substr(print_fasta, 1,nchar(str_extract(print_fasta, ">(.*)\\]")))
    print_fasta_seq = substr(print_fasta, 1+nchar(str_extract(print_fasta, ">(.*)\\]")),nchar(print_fasta))
  }
  
  cat(c(print_fasta_name,print_fasta_seq),
      sep = '\n')

} else{
  write(fasta_seq, paste0("P:/yonatan/saved.fasta.R/",fasta_name,".fasta"))
}  

}  
