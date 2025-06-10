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

#treatment = "mto1"
#test_df = "rnaseq"
#annotation="genes"
#context="CHG"

get_EC_4_KEGG <- function(test_df, annotation=NULL ,context=NULL) {
  ec2tair2product = read.csv("P:/TAIR10.1/EC2TAIR2product_index.csv")
  
  for (treatment in c("mto1","mto3")) {
    
    
    RNA_file.0 = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/met23/",treatment,"_vs_wt/all.transcripts.",treatment,"_vs_wt.DE.csv"))[,c("locus_tag","log2FoldChange","padj","EC.number","product")]
    RNA_file.0 = RNA_file.0[!(is.na(RNA_file.0$locus_tag) & is.na(RNA_file.0$EC.number) & is.na(RNA_file.0$product)), ]
    
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
      corr_file$padj = 0
      
      RNA_file.2 = RNA_file[,-grep("log2FoldChange|padj", names(RNA_file))]
      RNA_file_all = merge.data.frame(corr_file, RNA_file.2, by = "locus_tag")
      
      } else if (test_df == "meth") {
        if (context == "all") {
          meth_file = rbind(
            read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/genome_annotation/CG/",annotation,"_CG_genom_annotations.csv"))[,c("locus_tag","proportion1","proportion2","pValue")],
            read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/genome_annotation/CHG/",annotation,"_CHG_genom_annotations.csv"))[,c("locus_tag","proportion1","proportion2","pValue")],
            read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/genome_annotation/CHH/",annotation,"_CHH_genom_annotations.csv"))[,c("locus_tag","proportion1","proportion2","pValue")]
          )
        } else {
          meth_file = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/genome_annotation/",context,"/",annotation,"_",context,"_genom_annotations.csv"))[,c("locus_tag","proportion1","proportion2","pValue")]
        }
        
        meth_file$log2FoldChange = log2(meth_file$proportion2 / meth_file$proportion1)
        meth_file$padj = meth_file$pValue
        meth_file = meth_file[,-grep("proportion1|proportion2|pValue", names(meth_file))]
        
        
        
        RNA_file.2 = RNA_file[,-grep("log2FoldChange|padj", names(RNA_file))]
        RNA_file_all = merge.data.frame(meth_file, RNA_file.2, by = "locus_tag")
        
    } else if (test_df == "rnaseq") {
      RNA_file_all = RNA_file
      RNA_file_all$padj[is.na(RNA_file_all$padj)] = 0.999
    }
    
      
      RNA_loc_tmp = RNA_file_all[,grep("locus_tag|log2FoldChange|padj",names(RNA_file_all))]
      RNA_ec_tmp = RNA_file_all[,grep("EC.number|log2FoldChange|padj",names(RNA_file_all))]
      RNA_product_tmp = RNA_file_all[,grep("product|log2FoldChange|padj",names(RNA_file_all))]
      
      
      merged_tair = merge.data.frame(ec2tair2product, RNA_loc_tmp, by = "locus_tag")
      merged_ec = merge.data.frame(ec2tair2product, RNA_ec_tmp, by = "EC.number")
      merged_names = merge.data.frame(ec2tair2product, RNA_product_tmp, by = "product")
      
      df_2_save.0 = rbind(merged_tair[,c("EC.number","log2FoldChange","padj")],
                          merged_ec[,c("EC.number","log2FoldChange","padj")],
                          merged_names[,c("EC.number","log2FoldChange","padj")])
    
    
      df_2_save_sig = df_2_save.0[df_2_save.0$padj < 0.05,] %>% na.omit()
      
      MINlog2FC = ifelse(test_df == "corr", 0, ifelse(test_df == "meth", 0.1, 0.6))
      df_2_save_up = df_2_save_sig[df_2_save_sig$log2FoldChange >= MINlog2FC, c("EC.number","log2FoldChange")] %>% arrange(desc(log2FoldChange)) %>% distinct(EC.number, .keep_all = TRUE)
      df_2_save_down = df_2_save_sig[df_2_save_sig$log2FoldChange <= -MINlog2FC, c("EC.number","log2FoldChange")] %>% arrange(log2FoldChange) %>% distinct(EC.number, .keep_all = TRUE)
      
      df_2_save_both = merge.data.frame(df_2_save_up, df_2_save_down, by = "EC.number")
      df_2_save_both$log2FoldChange = paste(df_2_save_both$log2FoldChange.x,df_2_save_both$log2FoldChange.y, sep = " ")
      df_2_save_both = df_2_save_both[,-(2:3)]

      if (max(df_2_save.0$padj) >= 0.05) {
              df_2_save_all = df_2_save.0[df_2_save.0$padj >= 0.05, c("EC.number","log2FoldChange")] %>% na.omit()
      df_2_save_all$log2FoldChange = 0
      df_2_save = rbind(df_2_save_both,df_2_save_up,df_2_save_down,df_2_save_all) %>% distinct(EC.number, .keep_all = TRUE)
      } else {
        df_2_save = rbind(df_2_save_both,df_2_save_up,df_2_save_down) %>% distinct(EC.number, .keep_all = TRUE)
      }
      
      df_2_save$log2FoldChange[df_2_save$log2FoldChange == 0] = "#b0aea2"
      
      # normalize log2FC values to hex colors
      #normalized_values <- (df_2_save$log2FoldChange - min(df_2_save$log2FoldChange)) / (max(df_2_save$log2FoldChange) - min(df_2_save$log2FoldChange))
      #normalized_values <- (normalized_values * 2) - 1
      #gradient_func <- colorRampPalette(c("#5e66b8", "#e3dfc5", "#ad3232"))
      #df_2_save$log2FoldChange <- gradient_func(100)[as.integer((normalized_values + 1) * 49.5) + 1]

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
