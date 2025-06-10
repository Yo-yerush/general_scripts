library(KEGGREST)
library(dplyr)

for (treatment in c("mto1","mto3")) {
  
  tair2ec = read.csv("P:/TAIR10.1/uniProt_description_files/uniprotkb_arabidopsis_2023_11_13_short.txt",sep = "\t")[,c("TAIR","EC.number")]
  
  RNA_file = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/met23/",treatment,"_vs_wt/all.transcripts.",treatment,"_vs_wt.DE.csv"))[,c("log2FoldChange","padj","EC.number","product")]
  RNA_sig.0 = RNA_file[RNA_file$padj < 0.05,] %>% na.omit()
  
  # separate to a new raw by pattern
  RNA_sig.1 = as.data.frame(do.call(rbind, apply(RNA_sig.0, 1, function(x) {
    do.call(expand.grid, strsplit(x, ";"))
  })))
  RNA_sig.1$EC.number = as.character(RNA_sig.1$EC.number)
  
  for (i in 1:nrow(RNA_sig.1)) {
    if (length(grep("-",RNA_sig.1$EC.number[i])) != 0) {
      RNA_sig.1$EC.number[i] = gsub("- ","",paste(RNA_sig.1$EC.number[i],1:200, collapse = "; "))
    }
  }
  
  RNA_sig = as.data.frame(do.call(rbind, apply(RNA_sig.1, 1, function(x) {
    do.call(expand.grid, strsplit(x, ";"))
  })))
  RNA_sig$log2FoldChange = as.numeric(as.character(RNA_sig$log2FoldChange))
  RNA_sig$EC.number = as.character(RNA_sig$EC.number)
  
  RNA_up = RNA_sig[RNA_sig$log2FoldChange > 0.5,]
  RNA_down = RNA_sig[RNA_sig$log2FoldChange < -0.5,]
  RNA_both = merge.data.frame(RNA_up,RNA_down, by = "EC.number")
  
  ec_ids <- data.frame(EC.number = keggLink("enzyme", "map00270"))
  
 # ec_names = lapply(ec_ids$EC.number, function(ec) {
#    enzyme_details = keggGet(ec)
#    if (length(enzyme_details) > 0 && !is.null(enzyme_details[[1]]$NAME)) {
#      return(enzyme_details[[1]]$NAME)
#    } else {
#      return(NA)
#    }
#  })
#  
#  for (ii in 1:nrow(ec_ids)) {
#      ec_ids$product[ii] = paste(ec_names[[ii]], collapse = "")
#  }

  ec_ids$EC.number = gsub("ec:","",ec_ids$EC.number)
  
  ec_ids = as.data.frame(do.call(rbind, apply(ec_ids, 1, function(x) {
    do.call(expand.grid, strsplit(x, ";"))
  })))
  
  loop_list = list(
    list(RNA_both, "#f5958e #7a80c2"),#,"up&down",treatment))
    list(RNA_up, "#f5958e"),#,"up_regulated",treatment),
    list(RNA_down, "#7a80c2"),#,"down_regulated",treatment),
    list(RNA_file, "#d9deab")#,"non_sig",treatment),
  )
  
  df_2_save = data.frame(ec_number=NA,color=NA)
  iii=1
  
  for (iv in 1:4) {
    RNA = loop_list[[iv]][[1]]
    color = loop_list[[iv]][[2]]
    
    merged_ec = merge.data.frame(ec_ids,RNA, by = "EC.number")
    #merged_ec_names = merge.data.frame(ec_ids,RNA, by = "product")
    
    uni_ec = unique(merged_ec$EC.number)
    for (ec_number in uni_ec) {
      df_2_save[iii,] = c(ec_number,color)
      iii=iii+1
    }
  }
  
  df_2_save = df_2_save[!duplicated(df_2_save$ec_number),]
  write.table(df_2_save, paste0("P:/yonatan/methionine/rnaseq_23/EC_kegg_maps/",treatment,"_EC_for_kegg_maps.txt"), row.names = F, col.names = F, sep = " ", quote = F)
}

