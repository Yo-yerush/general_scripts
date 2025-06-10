align_fun = function(alignment.name,
                     protein_id_list,
                     FASTA_file = NULL, # FASTA file instead of IDs list
                     ordered_by_FASTA_file = F,
                     full_NCBI_name = T,
                     alignment.rows = NULL, # select row order
                     position, # position = c(start,end)
                     shading_Colors = "blues", #"blues", "reds", "greens", "grays", and "black" 
                     show_logo = F,
                     show_consensus = T, # when 'show_consensus = T', consensus will appear if there is *more than two seuence* to align
                     save_fasta_file = F,
                     TEX_only = F, # output only TEX file
                     insert_TEX_code = NA, # insert_TEX_code = c("code_1","code_2")
                     TEXshade_example_code = F, # print two latex codes for shading
                     TEX_open_script = F, # save TEX file and open script for edit
                     TEX_2_PDF = F, # convert TEX to PDF
                     path_to_save_pdf = "P:/yonatan/saved.pics.R/alignmemts/") {
  
  library(msa)
  
  if (TEXshade_example_code == T) {
    message("
######################################

# to shade all block (all rows)
\\shadeblock{tamplate_row}{start..end or MOTIF}{letter color}{shade color}

# to shade a selected row
\\shaderegion{tamplate_row}{start..end or MOTIF}{letter color}{shade color}


* tamplate_row >>> row number for the tamplate sequence
* start..end >>> position to shade by the 'tamplate_row' sequence
* MOTIF >>> AA sequence instead of position numbers
* letter/shade color >>> color name. first letter must be capital


examples:
\\shadeblock{1}{20..100}{White}{Black}
\\shadeblock{3}{ABCD}{Black}{Yellow}
\\shaderegion{2}{EFG}{White}{Green}

######################################
")
    stop("no file was saved", call. = F)
  }
  
  if (TEX_2_PDF == T) {  # after editing a TEX script only
    tools::texi2pdf(paste0(alignment.name,".tex"), clean=TRUE)
    stop()
  }

  original_path = getwd()
  setwd(path_to_save_pdf)
  
  alignment.name = gsub(" ","_", alignment.name)
  
#############################
  if (is.null(FASTA_file) == T) {
    aa = rentrez::entrez_fetch(db = "protein", id = protein_id_list, rettype = "fasta")
    #  gsub(">(.*)\\[Punica granatum\\]",">aa_alignment",aa,)
    write(aa,paste0(path_to_save_pdf,alignment.name,".fasta"))
    
    f_file = paste0(path_to_save_pdf,alignment.name,".fasta")
  } else {
    f_file = FASTA_file
  }
#############################
  
  mySequences <- Biostrings::readAAStringSet(f_file)
  
  if (full_NCBI_name == F) {
    names(mySequences) = substr(names(mySequences),1,14)
  }
  
  
  if (ordered_by_FASTA_file == F) {
    aa_alignment = msa(mySequences)
  } else {
    aa_alignment = msa(mySequences, order = "input")
    }
  
  
  if (show_consensus == T) {
    if (length(mySequences) > 2) {
      show_consensus = "bottom"
    } else {
      show_consensus = "none"
    }
  } else {show_consensus = "none"}

  
#  print(aa_alignment, show="complete")
  
  if (show_logo == T) {
    show_logo = "top"
  } else {
    show_logo = "none"
  }
  

    try({
      msa::msaPrettyPrint(aa_alignment, y = position,
                          output="pdf", showNames="left",
                          shadingColors = shading_Colors,
                          showLogo = show_logo, showConsensus = show_consensus,
                          verbose=F, askForOverwrite = T,
                          file = paste0(alignment.name,".pdf"),
                          subset = alignment.rows,
                          furtherCode = insert_TEX_code, consensusColors = "Gray")
    })

  
  gc()
  
#  if (save_fasta_file == F) {
#    file.remove(paste0(path_to_save_pdf,alignment.name,".fasta"))
#  }
  
  if (TEX_open_script == T) {
    file.edit(paste0("P:/yonatan/saved.pics.R/alignmemts/",alignment.name,".tex"))
    stop()
  } 
  
  if (TEX_only == F) {
    try({
    tools::texi2pdf(paste0(alignment.name,".tex"), clean=TRUE)
    file.remove(paste0(alignment.name,".tex"))
    stop()
    })
  }

  setwd(original_path)
}
