rename_seq_TAIR.2.TAIR <- function(gr_obj, is_TAIR10.1) {
  #### change seqnames function
  gr_obj_trimmed = gr_obj
  
  if (is_TAIR10.1) {
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("Chr1","NC_003070.9",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("Chr2","NC_003071.7",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("Chr3","NC_003074.8",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("Chr4","NC_003075.7",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("Chr5","NC_003076.8",seqlevels(gr_obj_trimmed)))
  } else {
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003070.9","Chr1",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003071.7","Chr2",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003074.8","Chr3",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003075.7","Chr4",seqlevels(gr_obj_trimmed)))
    gr_obj_trimmed = renameSeqlevels(gr_obj_trimmed, gsub("NC_003076.8","Chr5",seqlevels(gr_obj_trimmed)))
  }
  
  return(gr_obj_trimmed)
}