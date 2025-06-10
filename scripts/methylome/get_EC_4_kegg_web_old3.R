library(KEGGREST)
library(dplyr)

if (F) {
  ec_ids <- data.frame(EC.number = keggLink("enzyme", "map00270"))
  # call gene names by EC --- take some time...
  product2ec = lapply(ec_ids$EC.number, function(ec) {
    enzyme_details = keggGet(ec)
    if (length(enzyme_details) > 0 && !is.null(enzyme_details[[1]]$NAME)) {
      return(enzyme_details[[1]]$NAME)
    } else {
      return(NA)
    }
  })
  for (ii in 1:nrow(ec_ids)) {
    ec_ids$product[ii] = paste(product2ec[[ii]], collapse = "")
  }
  tair_ids <- keggLink("ath", ec_ids$EC.number)
  ec2tair = data.frame(locus_tag = tair_ids, EC.number = names(tair_ids))
  ec2tair2product.0 = merge.data.frame(ec2tair, ec_ids, by = "EC.number", all.y = T)
  
  ec2tair2product.0$locus_tag = gsub("ath:","",ec2tair2product.0$locus_tag)
  ec2tair2product.0$EC.number = gsub("ec:","",ec2tair2product.0$EC.number)
  
  ec2tair2product = as.data.frame(do.call(rbind, apply(ec2tair2product.0, 1, function(x) {
    do.call(expand.grid, strsplit(x, ";"))
  })))
  
  ec2tair2product$locus_tag = as.character(ec2tair2product$locus_tag)
  ec2tair2product$EC.number = as.character(ec2tair2product$EC.number)
  ec2tair2product$product = tolower(as.character(ec2tair2product$product))
  
  write.csv(ec2tair2product, "P:/TAIR10.1/EC2TAIR2product_index.csv", row.names = F)
}

get_EC_4_KEGG <- function(test_df, annotation=NULL ,context=NULL) {
  ec2tair2product = read.csv("P:/TAIR10.1/EC2TAIR2product_index.csv")
  
  for (treatment in c("mto1","mto3")) {
    
    
    RNA_file.0 = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/met23/",treatment,"_vs_wt/all.transcripts.",treatment,"_vs_wt.DE.csv"))[,c("locus_tag","log2FoldChange","padj","EC.number","product")]
    RNA_file.0$product = tolower(RNA_file.0$product)
    RNA_file.0$product = gsub(" [0-9]+$", "", RNA_file.0$product)
    RNA_file.0$product = gsub(" \\(.*?\\)", "", RNA_file.0$product)
    RNA_file.0$product = gsub(" family protein", "", RNA_file.0$product)
    RNA_file.0$product = gsub(" family protein", "", RNA_file.0$product)
    RNA_file.0$product = gsub("-like protein", "", RNA_file.0$product)
    RNA_file.0$product = gsub("-like", "", RNA_file.0$product)
    
    RNA_file.with_ec = RNA_file.0[grep("[0-9+]",RNA_file.0$EC.number),]
    RNA_file.without_ec = RNA_file.0[-grep("[0-9+]",RNA_file.0$EC.number),]
    
    RNA_file.with_ec = as.data.frame(do.call(rbind, apply(RNA_file.with_ec, 1, function(x) {
      do.call(expand.grid, strsplit(x, "; "))
    })))
    RNA_file.with_ec$locus_tag = as.character(RNA_file.with_ec$locus_tag)
    RNA_file.with_ec$log2FoldChange = as.numeric(as.character(RNA_file.with_ec$log2FoldChange))
    RNA_file.with_ec$padj = as.numeric(as.character(RNA_file.with_ec$padj))
    RNA_file.with_ec$EC.number = as.character(RNA_file.with_ec$EC.number)
    RNA_file.with_ec$product = as.character(RNA_file.with_ec$product)
    
    RNA_file = rbind(RNA_file.with_ec,RNA_file.without_ec)
    
    if (test_df == "corr") {
      if (context == "all") {
        corr_file = rbind(
          read.csv(paste0("P:/yonatan/methionine/rnaseq_23/corr_with_methylations/",treatment,"/CG/",annotation,".corr.CG.",treatment,".csv"))[,c("locus_tag","cor","pval")],
          read.csv(paste0("P:/yonatan/methionine/rnaseq_23/corr_with_methylations/",treatment,"/CHG/",annotation,".corr.CHG.",treatment,".csv"))[,c("locus_tag","cor","pval")],
          read.csv(paste0("P:/yonatan/methionine/rnaseq_23/corr_with_methylations/",treatment,"/CHH/",annotation,".corr.CHH.",treatment,".csv"))[,c("locus_tag","cor","pval")]
        )
      } else {
        corr_file = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/corr_with_methylations/",treatment,"/",context,"/",annotation,".corr.",context,".",treatment,".csv"))[,c("locus_tag","cor","pval")]
      }
  
      names(corr_file)[2:3] = c("log2FoldChange", "padj")
      
      RNA_file.2 = RNA_file[,-grep("log2FoldChange|padj", names(RNA_file))]
      RNA_file_all = merge.data.frame(corr_file, RNA_file.2, by = "locus_tag")
      RNA_sig = RNA_file_all
      
      
      } else if (test_df == "meth") {
        if (context == "all") {
          meth_file = rbind(
            read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/genome_annotation/CG/",annotation,"_CG_genom_annotations.csv"))[,c("locus_tag","regionType")],
            read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/genome_annotation/CHG/",annotation,"_CHG_genom_annotations.csv"))[,c("locus_tag","regionType")],
            read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/genome_annotation/CHH/",annotation,"_CHH_genom_annotations.csv"))[,c("locus_tag","regionType")]
          )
        } else {
          meth_file = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/genome_annotation/",context,"/",annotation,"_",context,"_genom_annotations.csv"))[,c("locus_tag","regionType")]
        }
        
        names(meth_file)[2] = c("log2FoldChange")
        meth_file$log2FoldChange[grep("gain", meth_file$log2FoldChange)] = 1
        meth_file$log2FoldChange[grep("loss", meth_file$log2FoldChange)] = -1
        meth_file$log2FoldChange = as.numeric(meth_file$log2FoldChange)
        
        RNA_file.2 = RNA_file[,-grep("log2FoldChange", names(RNA_file))]
        RNA_file_all = merge.data.frame(meth_file, RNA_file.2, by = "locus_tag")
        RNA_sig = RNA_file_all
        
        
    } else if (test_df == "rnaseq") {
      RNA_sig = RNA_file[RNA_file$padj < 0.05,] %>% na.omit()
      RNA_file_all = RNA_file
    }
    
    RNA_up = RNA_sig[RNA_sig$log2FoldChange >= 0.6,]
    RNA_down = RNA_sig[RNA_sig$log2FoldChange <= -0.6,]
    #RNA_both = merge.data.frame(RNA_up,RNA_down, by = "locus_tag")
    
    
    loop_list = list(
      #list(RNA_both, "#f5958e #7a80c2"),#,"up&down",treatment))
      list(RNA_up, "#f5958e"),#,"up_regulated",treatment),
      list(RNA_down, "#7a80c2"),#,"down_regulated",treatment),
      list(RNA_file_all, "#d9deab")#,"non_sig",treatment),
    )
    
    df_2_save.0 = data.frame(EC.number=NA, log2FoldChange=NA)
    
    for (iv in 1:3) {
      RNA = loop_list[[iv]][[1]]
      color = loop_list[[iv]][[2]]
      
      #RNA_loc_tmp = data.frame(locus_tag = paste0(unlist(RNA[,grep("locus_tag",names(RNA))]))) # why unlist? remove it and add 'log2FC' column
      #RNA_ec_tmp = data.frame(EC.number = paste0(unlist(RNA[,grep("EC.number",names(RNA))])))
      #RNA_product_tmp = data.frame(product = paste0(unlist(RNA[,grep("product",names(RNA))])))
      
      RNA_loc_tmp = RNA[,grep("locus_tag|log2FoldChange",names(RNA))]
      RNA_ec_tmp = RNA[,grep("EC.number|log2FoldChange",names(RNA))]
      RNA_product_tmp = RNA[,grep("product|log2FoldChange",names(RNA))]
      
      
      merged_tair = merge.data.frame(ec2tair2product, RNA_loc_tmp, by = "locus_tag")
      merged_ec = merge.data.frame(ec2tair2product, RNA_ec_tmp, by = "EC.number")
      merged_names = merge.data.frame(ec2tair2product, RNA_product_tmp, by = "product")
      
      df_2_save.0 = rbind(df_2_save.0,
                          merged_tair[,c("EC.number","log2FoldChange")],
                          merged_ec[,c("EC.number","log2FoldChange")],
                          merged_names[,c("EC.number","log2FoldChange")])

    }
    
    df_2_save_up = df_2_save.0[df_2_save.0$color == "#f5958e",]
    df_2_save_down = df_2_save.0[df_2_save.0$color == "#7a80c2",]
    df_2_save_both = merge.data.frame(df_2_save_up, df_2_save_down, by = "ec_number")
    df_2_save_both$color = paste(df_2_save_both$color.x,df_2_save_both$color.y, sep = " ")
    df_2_save_both = df_2_save_both[,-c(2:3)]
    df_2_save_all = df_2_save.0[df_2_save.0$color == "#d9deab",]
    
    df_2_save = rbind(df_2_save_both,df_2_save_up,df_2_save_down,df_2_save_all)
    df_2_save = df_2_save[!duplicated(df_2_save$ec_number),]
    
    if (test_df == "corr") {
      file_suffix = paste0("_corr_",annotation,"_",context,"_EC_for_kegg_maps.txt")
    } else if (test_df == "meth") {
      file_suffix = paste0("_meth_",annotation,"_",context,"_EC_for_kegg_maps.txt")
    } else if (test_df == "rnaseq") {
      file_suffix = "_rnaseq_EC_for_kegg_maps.txt"
    }

    write.table(df_2_save, paste0("P:/yonatan/methionine/rnaseq_23/EC_kegg_maps/",treatment,file_suffix), row.names = F, col.names = F, sep = " ", quote = F)
  }
}

get_EC_4_KEGG("rnaseq")

for (annotation.l in c("Genes","Promoters")) {
  for (context.l in c("CG","CHG","CHH","all")) {
    get_EC_4_KEGG("meth", annotation.l, context.l)
  }
}

for (annotation.ll in c("genes","promoters")) {
  for (context.ll in c("CG","CHG","CHH","all")) {
    get_EC_4_KEGG("corr", annotation.ll, context.ll)
  }
}