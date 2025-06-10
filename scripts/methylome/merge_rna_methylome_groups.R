library(dplyr)
#library(tidyr)
library(stringr)
library(purrr)
library(writexl)
library(openxlsx)

#make_it_unique = F # if TRUE, each gene can obtain in one group only (the first by order)

offspring_fun <- function(go_id, xx = as.list(GO.db::GOBPOFFSPRING)) { # 'GOBPCHILDREN' for child terms
  
  child_terms_0 = as.character(xx[[go_id]])
  child_terms = child_terms_0
  
  for (i in 1:length(child_terms_0)) {
    child_terms = c(child_terms, as.character(xx[[child_terms[i]]]))
  } 
  
  return(child_terms[!is.na(child_terms)] %>% unique()) # %>% paste(collapse = "|"))
}


################

grep_position <- function(x) {
  vec = NULL
  for (terms_l in x) {
    vec = c(vec, grep(terms_l, all_res$GO.biological.process))
  }
  return(unique(vec))
}


for (treatment in c("mto1_vs_wt", "mto3_vs_wt", "dCGS_vs_EV", "SSE_H_vs_EV", "SSE_L_vs_EV")) {
  for (make_it_unique in c(F)) {
    unique_or_not = ifelse(make_it_unique, "_unique","")
    
    all_res = as.data.frame(readxl::read_xlsx("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_genes.xlsx", sheet = treatment))
    #all_res = all_res[,-ncol(all_res)] # %>% relocate("transcript_id", .after = "gene_id")
    
    ################
    # RdDM pathway
    rddm = rbind(read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_RNA-directed DNA methylation_an epigenetic pathway of increasing complexity/tairs_ID.txt",
                            sep = "\t", header = T),
                 read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_Non-canonical RNA-directed DNA methylation/tairs_ID.txt",
                            sep = "\t", header = T)
    ) %>% distinct()
    
    names(rddm) = "gene_id"
    
    rddm_mto1 = merge.data.frame(rddm, all_res, by = "gene_id") %>% arrange(RNA_padj)
    ################
    
    ################
    # Histone Lysine Methyltransferases
    HLM = rbind(read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_Histone Lysine Methyltransferases/tairs_ID.txt",
                           sep = "\t", header = T),
                read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_Plant SET Domain-containing Proteins_Structure, Function and Regulation/tairs_id.txt",
                           sep = "\t", header = T)) %>% 
      distinct()
    
    names(HLM)[1] = "gene_id"
    HLM_mto1 = merge.data.frame(HLM, all_res, by = "gene_id") %>%
      arrange(RNA_padj)
    ################
    
    ################
    # Royal Family Proteins
    RF = read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/Royal Family proteins/At_Agenet_Tudor_family.txt",
                    sep = "\t", header = T)
    names(RF)[1] = "gene_id"
    RF$gene_id = gsub("DUF\\d+","", RF$gene_id)
    RF_mto1 = merge.data.frame(RF, all_res, by = "gene_id") %>% dplyr::select(-Agenet.Tudor.Domain, -Other.Domain)
    ################
    
    ################
    # Cohen SSE
    sse = read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/Cohen_14 SSE/gene_list.txt",
                     sep = "\t", header = T)
    names(sse)[1] = "gene_id"
    sse$gene_id = toupper(sse$gene_id)
    sse_mto1 = merge.data.frame(sse, all_res, by = "gene_id") %>% dplyr::select(-Title)
    ################
    # GO term offsprings
    child_terms_epigenetic = offspring_fun("GO:0040029") # epigenetic regulation of_gene expression
    child_terms_chromatin_org = offspring_fun("GO:0006325") # chromatin organization
    child_terms_chromatin_rem = offspring_fun("GO:0006338") # chromatin remodeling
    
    ################
    xl_0_list = list(
      DNA_methyltransferase = rbind(all_res[grep("^2\\.1\\.1\\.37", all_res$EC),], # EC_2.1.1.37
                                    all_res[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(all_res$Protein.names)),], # DNA_C5_MT
                                    all_res[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(all_res$short_description)),], # DNA_C5_MT
                                    all_res[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(all_res$Computational_description)),],
                                    all_res[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(all_res$GO.biological.process)),]# DNA_C5_MT
      ),

      Histone_Lysine_MTs = rbind(HLM_mto1, # %>% filter(RNA_pvalue < 0.05),
                                 all_res[grep("SET domain",  gsub("SET-domain","SET domain",all_res$short_description)),], # %>% filter(RNA_pvalue < 0.05), # SET_domain
                                 #all_res[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(all_res$Gene_description)),], # %>% filter(RNA_pvalue < 0.05), # SDG
                                 all_res[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(all_res$Symbol)),], # %>% filter(RNA_pvalue < 0.05), # SDG
                                 all_res[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(all_res$old_symbols)),], # %>% filter(RNA_pvalue < 0.05), # SDG
                                 all_res[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(all_res$Protein.names)),], # %>% filter(RNA_pvalue < 0.05), # SDG
                                 all_res[grep("class v-like sam-binding methyltransferase", tolower(all_res$Protein.families)),] # %>% filter(RNA_pvalue < 0.05) # SAM_MT
      ),
      
      RdDM_pathway = rddm_mto1, # %>% filter(RNA_pvalue < 0.05),
      Royal_Family_Proteins = RF_mto1, # %>% filter(RNA_pvalue < 0.05),
      
      #DNA_methylation_related = rbind(all_res[grep("dna methylation", tolower(all_res$Function..CC.)),] %>% filter(RNA_pvalue < 0.05), # DNA_methylation_function
      #                                all_res[grep("dna methylation", tolower(all_res$Gene.Ontology..biological.process.)),] %>% filter(RNA_pvalue < 0.05), # DNA_methylation_bp
      #                                all_res[grep("rddm", tolower(all_res$Function..CC.)),] %>% filter(RNA_pvalue < 0.05) # RdDM_function
      #),
      chromatin_remodeling = rbind(all_res[grep("chromatin remodeling|chromatin remodeler", tolower(all_res$GO.biological.process)),], # %>% filter(RNA_pvalue < 0.05), # chromatin_remodeling_BP
                                   all_res[grep("chromatin remodeling|chromatin remodeler", tolower(all_res$Gene_description)),], # %>% filter(RNA_pvalue < 0.05), # chromatin_remodeling_pro_name
                                   all_res[grep("chromatin remodeling|chromatin remodeler", tolower(all_res$Computational_description)),], # %>% filter(RNA_pvalue < 0.05), # chromatin_remodeling_pro_name
                                   all_res[grep("chromatin remodeling|chromatin remodeler", tolower(all_res$Function)),], # %>% filter(RNA_pvalue < 0.05), # chromatin_remodeling_pro_name
                                   all_res[grep("chromatin remodeling|chromatin remodeler", tolower(all_res$short_description)),] # %>% filter(RNA_pvalue < 0.05) # chromatin remodeling_function
      ),
      #other_methylation = rbind(all_res[grep("class i-like sam-binding methyltransferase", tolower(all_res$Protein.families)),], # %>% filter(RNA_pvalue < 0.05), # SAM_MT
      #                          all_res[grep("class iv-like sam-binding methyltransferase", tolower(all_res$Protein.families)),] # %>% filter(RNA_pvalue < 0.05)#, # SAM_MT
      #                          #all_res[grep("^2\\.1\\.1\\.", all_res$EC.number),] %>% filter(RNA_pvalue < 0.05) # EC_2.1.1
      #),
      chromatin_organization_related = rbind(all_res[grep("chromosome", tolower(all_res$GO.cellular.component)),]  %>% filter(RNA_pvalue < 0.05), # chromosome_CC
                                      all_res[grep("chromatin", tolower(all_res$GO.cellular.component)),]  %>% filter(RNA_pvalue < 0.05), # chromatin_CC
                                      all_res[grep("chromatin", tolower(all_res$Function)),] %>% filter(RNA_pvalue < 0.05), # chromatin_function
                                      all_res[grep("nucleosome", tolower(all_res$GO.cellular.component)),]  %>% filter(RNA_pvalue < 0.05), # nucleosome_CC
                                      #all_res[grep("zinc finger", tolower(all_res$Protein.names)),],  %>% filter(RNA_pvalue < 0.05), # zinc_finger_pro_name
                                      all_res[grep("jumonji", tolower(all_res$short_description)),]  %>% filter(RNA_pvalue < 0.05), # jumonji
                                      all_res[grep("helicase domain", tolower(all_res$short_description)),]  %>% filter(RNA_pvalue < 0.05), # helicase_domain
                                      all_res[grep("histone h3 acetylation",tolower(all_res$GO.biological.process)),] %>% filter(RNA_pvalue < 0.05), # histone_H3_acetylation_BP
                                      all_res[grep_position(child_terms_chromatin_org),]  %>% filter(RNA_pvalue < 0.05),
                                      all_res[grep_position(child_terms_chromatin_rem),]  %>% filter(RNA_pvalue < 0.05)
      ),
      
      #Cohen_SSE = sse_mto1, # %>% filter(RNA_pvalue < 0.05),
      
      # keep al rows the contain at least one DMR in list of genes (from GO ID)
      epigenetic_reg._gene_expression = all_res[grep_position(child_terms_epigenetic),] # %>% 
        #filter(!(is.na(CG_Genes) & is.na(CHG_Genes) & is.na(CHH_Genes) & is.na(CG_Promoters) & is.na(CHG_Promoters) & is.na(CHH_Promoters)))
    )
    
    
    # empty df with correct headers
    tmp_df = data.frame(matrix(ncol = length(names(xl_0_list[[1]])), nrow = 0))
    names(tmp_df) = names(xl_0_list[[1]])
    
    # create data frame to make it unique by rows
    for (i in 1:length(xl_0_list)) {
      tmp_loop_df = xl_0_list[[i]]
      tmp_loop_df$tmp1 = names(xl_0_list)[i]
      
      ### make it unique or not (if each gene will be in one group or all related groups)
      tmp_loop_df$tmp2 = paste0(tmp_loop_df$gene_id,tmp_loop_df$CG_Genes,tmp_loop_df$CHG_Genes,tmp_loop_df$CHH_Genes,tmp_loop_df$CG_Promoters,tmp_loop_df$CHG_Promoters,tmp_loop_df$CHH_Promoters)
      if (make_it_unique == F) {
        tmp_loop_df = tmp_loop_df %>% distinct(tmp2, .keep_all = TRUE) %>% dplyr::select(-tmp2)
      }
      
      tmp_df = rbind(tmp_df, tmp_loop_df)
    }
    
    if (make_it_unique) {
      tmp_df = tmp_df %>% distinct(tmp2, .keep_all = TRUE) %>% dplyr::select(-tmp2)
    }
    
    # split again to a list
    xl_list <- split(tmp_df, tmp_df$tmp1)
    xl_list = xl_list[unique(tmp_df$tmp1)] # order it as suppose to be...
    
    ################################################################################
    
    ################################################################################
    remove_dup_DMR <- function(y) {
      y <- as.character(unique(unlist(strsplit(y, ","))))
      paste(y, collapse = ",")
    }
    
    
    ################ eddit xl_list data frames
    xl_list <- lapply(xl_list, function(x) {
      
      # Group by gene_id and summarize (add DMRs values for the gene, make in onw row for each gene)
      x = x %>%
        group_by(gene_id) %>%
        summarise(
          across(contains("_Genes") | contains("_Promoters"), ~remove_dup_DMR(paste(., collapse = ","))),  # apply remove_dup_DMR function for DMR columns
          across(!contains("_Genes") & !contains("_Promoters"), first)  # for other columns
        ) %>%
        as.data.frame() %>%
        mutate(across(contains("_Genes") | contains("_Promoters"), ~ str_replace_all(.x, "NA", ""))) %>%
        mutate(across(contains("_Genes") | contains("_Promoters"), ~ str_replace_all(.x, ",", ", "))) %>% # change delimiter
        #mutate(gene = replace_na(gene, "")) %>%
        
        relocate(contains("_Genes") | contains("_Promoters"), .before = Symbol) %>%
        #relocate("Function..CC.", .after = "short_description") %>%
        #relocate("gene_model_type", .before = "DOI.ID") %>%
        select(-tmp1) %>%
        
        arrange(RNA_padj) 
      
      #######

      return(x)
    })
    
    ################################################################################
    
    # make only 'epigenetic' sheet unique with others
    dup_genes = lapply(xl_list[-length(xl_list)], function(x) {x[[1]]}) %>% unlist() %>% as.character() %>% unique()
      
    tmp_arg = xl_list[["epigenetic_reg._gene_expression"]]
    tmp_arg = tmp_arg[!tmp_arg$gene_id %in% dup_genes, ]
    xl_list[["epigenetic_reg._gene_expression"]] = tmp_arg
    
    ################################################################################
    
    ################ # add beck original columns of 'sse' and 'royal family proteins'
    #xl_list[["Cohen_SSE"]] = merge.data.frame(sse, xl_list[["Cohen_SSE"]], by = "gene_id") %>% 
    #  relocate("Title", .after = "gene") %>%
    #  arrange(RNA_padj) %>%
    #  arrange(Title)
    
    xl_list[["Royal_Family_Proteins"]] = merge.data.frame(RF, xl_list[["Royal_Family_Proteins"]], by = "gene_id") %>% 
      relocate(c("Agenet.Tudor.Domain", "Other.Domain"), .after = "Symbol") %>%
      arrange(RNA_padj)
    
    ################################################################################
    
    ################################################################################
    clean_ASCII <- function(x) {
      x = gsub("\002", " ", x)
      x = gsub("\036", " ", x)
      
      #  x = gsub("[[:punct:]]", " ", x)
      #x = iconv(x, from = 'UTF-8', to = 'ASCII')
      return(x)
    }
    ################
    xl_headers = names(xl_list[[1]])
    numeric_cols = grep("RNA_|CG_|CHG_|CHH_", xl_headers)
    RNA_cols = grep("RNA_", xl_headers)
    DMRs_cols = grep("CG_|CHG_|CHH_", xl_headers)
    p_cols = grep("RNA_p", xl_headers)
    lfc_cols = grep("RNA_log2FC", xl_headers)
    other_cols = (grep("CHH_Promoters", xl_headers)+2) : length(xl_headers)
    ################
    # save and edit EXCEL
    wb <- createWorkbook()
    # Define styles
    style_up <- createStyle(bgFill = "#f59d98")
    style_down <- createStyle(bgFill = "#c3ccf7")
    style_p <- createStyle(bgFill = "#f7deb0")
    style_other <- createStyle(bgFill = "#daf7d7")
    cell_n_font_style <- createStyle(border = "TopBottomLeftRight", borderColour = "black")
    header_style <- createStyle(textDecoration = "bold", border = "TopBottomLeftRight", borderStyle = "double")
    
    DMR_style_up <- createStyle(fgFill = "#f59d96", border = "TopBottomLeftRight", borderColour = "black")
    DMR_style_down <- createStyle(fgFill = "#c3ccf7", border = "TopBottomLeftRight", borderColour = "black")
    DMR_style_shared <- createStyle(fgFill = "#daf7d7", border = "TopBottomLeftRight", borderColour = "black")
    
    for (sheet_name in names(xl_list)) {
      
      df <- xl_list[[sheet_name]]
      df <- data.frame(lapply(df, clean_ASCII))
      
      #### for make this columns as numeric
      #duplicate_DMRs_values_rows = df[,DMRs_cols] %>% # find duplicate DMRs values (they cant be numeric...)
      #  mutate(across(everything(), ~grepl(",", .))) %>%
      #  reduce(`|`) %>%
      #  which()
      
      #non_duplicate_DMRs_values_rows = (1:nrow(df))[-duplicate_DMRs_values_rows]
      #df[-duplicate_DMRs_values_rows,DMRs_cols] = sapply(df[-duplicate_DMRs_values_rows, DMRs_cols], as.numeric)
      
      df[,RNA_cols] = sapply(df[,RNA_cols], as.numeric)
      #### 
      
      addWorksheet(wb, sheet_name)
      #  setColWidths(wb, sheet_name, cols = other_cols, widths = 10)
      #  setColWidths(wb, sheet_name, cols = lfc_cols, widths = 4)
      #  setColWidths(wb, sheet_name, cols = p_cols, widths = 6)
      
      
      writeData(wb, sheet_name, df)
      
      addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
      addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
      
      # colors for DMRs values (minus as blue as plus as red)
      for (i.c in DMRs_cols) {
        for (i.r in 1:nrow(df)) {
          
          # for cells with both plus and minus direction
          if ( length(grep("^-.*, [0-9]", df[i.r,i.c])) != 0 | length(grep("^[0-9].*, -[0-9]", df[i.r,i.c])) != 0 ) {
            addStyle(wb, sheet_name, style = DMR_style_shared, 
                     rows = i.r + 1, # first row in excel is the headers
                     cols = i.c,
                     gridExpand = FALSE)
            
            # for cells with minus direction
          } else if ( length(grep("-", df[i.r,i.c])) != 0 ) {
            addStyle(wb, sheet_name, style = DMR_style_down, 
                     rows = i.r + 1, # first row in excel is the headers
                     cols = i.c,
                     gridExpand = FALSE)
            
            # for cells with plus direction 
          } else if ( length(grep("^[0-9]", df[i.r,i.c])) != 0 ) {
            addStyle(wb, sheet_name, style = DMR_style_up, 
                     rows = i.r + 1, # first row in excel is the headers
                     cols = i.c,
                     gridExpand = FALSE)
          }
        }     
      }
      #      addStyle(wb, sheet_name, style = header_style, rows = duplicate_DMRs_values_rows , cols = DMRs_cols, gridExpand = TRUE)
      #      addStyle(wb, sheet_name, style = header_style, rows = non_duplicate_DMRs_values_rows , cols = DMRs_cols, gridExpand = TRUE)
      
      conditionalFormatting(wb, sheet_name, cols = p_cols[1], rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
      conditionalFormatting(wb, sheet_name, cols = p_cols[2], rows = 2:(nrow(df)+1), rule = "<0.05", style = style_p)
      
      for (col in lfc_cols) {
        conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = ">0", style = style_up)
        conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df)+1), rule = "<0", style = style_down)
      }
      
      for (col in other_cols[other_cols %% 2 == 0]) {
        conditionalFormatting(wb, sheet_name, style = style_other, rule = "!=0", rows = 2:(nrow(df)+1), cols = col, gridExpand = TRUE)
        conditionalFormatting(wb, sheet_name, style = style_other, rule = "==0", rows = 2:(nrow(df)+1), cols = col, gridExpand = TRUE)
      }
    }
    saveWorkbook(wb, paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/", treatment,unique_or_not,"_groups.xlsx"), overwrite = T)
    
  }
}
