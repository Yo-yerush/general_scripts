pro.id_2_fasta <- function(protein_id) {
  
  if (protein_id == "XP_031373666.1") {
    sdh_name = paste0("SDH_1")
  }
  if (protein_id == "XP_031373667.1") {
    sdh_name = paste0("SDH_1")
  }
  if (protein_id == "XP_031373668.1") {
    sdh_name = paste0("SDH_1")
  }
  if (protein_id == "XP_031373669.1") {
    sdh_name = paste0("SDH_1")
  }
  
  if (protein_id == "XP_031407288.1") {
    sdh_name = paste0("SDH_3.1")
  }
  if (protein_id == "XP_031377912.1") {
    sdh_name = paste0("SDH_3.2")
  }
  if (protein_id == "XP_031405648.1") {
    sdh_name = paste0("SDH_3a.1")
  }
  if (protein_id == "XP_031405900.1") {
    sdh_name = paste0("SDH_3a.2")
  }
  if (protein_id == "XP_031405901.1") {
    sdh_name = paste0("SDH_3a.2")
  }
  if (protein_id == "XP_031405902.1") {
    sdh_name = paste0("SDH_3a.2")
  }
  
  if (protein_id == "XP_031377057.1") {
    sdh_name = paste0("SDH_4")
  }
  
  
  library(msa)
  library(ape)
  library(stringr)
  
  aa = rentrez::entrez_fetch(db = "protein", id = protein_id, 
                             rettype = "fasta")
  aa.2 = substr(aa, 1+nchar(str_extract(aa, ">X(.*)\\]")),nchar(aa))
  
  aa.3 = paste0(">",sdh_name," (",protein_id,")",aa.2)
  print(aa.3)
}  
