library(dplyr)
# library(tidyr)
library(stringr)
library(purrr)
library(writexl)
library(openxlsx)
library(ggplot2)
library(cowplot)

# make_it_unique = F # if TRUE, each gene can obtain in one group only (the first by order)

output_res <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/genes_group_results/"


offspring_fun <- function(go_id, xx = as.list(GO.db::GOBPOFFSPRING)) { # 'GOBPCHILDREN' for child terms

  child_terms_0 <- as.character(xx[[go_id]])
  child_terms <- child_terms_0

  for (i in 1:length(child_terms_0)) {
    child_terms <- c(child_terms, as.character(xx[[child_terms[i]]]))
  }

  return(child_terms[!is.na(child_terms)] %>% unique()) # %>% paste(collapse = "|"))
}


################

grep_position <- function(x) {
  vec <- NULL
  for (terms_l in x) {
    vec <- c(vec, grep(terms_l, all_res$GO.biological.process))
  }
  return(unique(vec))
}


for (treatment in c("mto1_vs_wt", "mto3_vs_wt", "dCGS_vs_EV", "SSE_high_vs_EV", "SSE_low_vs_EV", "SSE_high_vs_SSE_low")) {
  for (make_it_unique in c(F)) {
    unique_or_not <- ifelse(make_it_unique, "_unique", "")

    all_res <- as.data.frame(readxl::read_xlsx("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/merged_results_mtos_all_genes.xlsx", sheet = treatment, progress = FALSE))
    # all_res = all_res[,-ncol(all_res)] # %>% relocate("transcript_id", .after = "gene_id"

    cat(paste0(treatment, "..."))
    ################
    # RdDM pathway
    rddm <- rbind(
      read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_RNA-directed DNA methylation_an epigenetic pathway of increasing complexity/tairs_ID.txt",
        sep = "\t", header = T
      ),
      read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_Non-canonical RNA-directed DNA methylation/tairs_ID.txt",
        sep = "\t", header = T
      )
    ) %>% distinct()

    names(rddm) <- "gene_id"

    rddm_mto1 <- merge.data.frame(rddm, all_res, by = "gene_id") %>% arrange(RNA_padj)
    ################

    ################
    # Histone Lysine Methyltransferases
    HLM <- rbind(
      read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_Histone Lysine Methyltransferases/tairs_ID.txt",
        sep = "\t", header = T
      ),
      read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_Plant SET Domain-containing Proteins_Structure, Function and Regulation/tairs_id.txt",
        sep = "\t", header = T
      )
    ) %>%
      distinct()

    names(HLM)[1] <- "gene_id"
    HLM_mto1 <- merge.data.frame(HLM, all_res, by = "gene_id") %>%
      arrange(RNA_padj)
    ################
    cat(".")
    ################
    # Royal Family Proteins
    RF <- read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/Royal Family proteins/At_Agenet_Tudor_family.txt",
      sep = "\t", header = T
    )
    names(RF)[1] <- "gene_id"
    RF$gene_id <- gsub("DUF\\d+", "", RF$gene_id)
    RF_mto1 <- merge.data.frame(RF, all_res, by = "gene_id") %>% dplyr::select(-Agenet.Tudor.Domain, -Other.Domain)
  
    ################
    # DNA de-Methylases
    DNA_deMTs <- all_res[grep("AT4G34060|AT5G04560|AT2G36490|AT3G10010", all_res$gene_id), ] %>%
      distinct(gene_id, .keep_all = T) %>%
      arrange(RNA_padj)

    # histone de-Methylases
    histone_deMTs <- all_res[grep("LDL|FLD|ELF|IBM|JMJ|REF", all_res$Symbol), ] %>%
      distinct(gene_id, .keep_all = T) %>%
      arrange(RNA_padj)
    ################
    # Ash leaf
    Ash_leaf <- read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_Girija_23 methylation Plant Physiol/Ash_leaves.txt",
      sep = "\t", header = T
    )

    Ash_merged <- merge.data.frame(Ash_leaf, all_res, by = "gene_id")

    ################
    # primary/secondary metabolism (https://doi.org/10.1093/gbe/evv217)
    primary_metabolism_0 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_Evolutionary Rate Heterogeneity of Primary and Secondary Metabolic Pathway Genes in Arabidopsis thaliana/primary_metabolism_genes_pathways_with_groups.csv")

    secondary_metabolism_0 <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/sup_Evolutionary Rate Heterogeneity of Primary and Secondary Metabolic Pathway Genes in Arabidopsis thaliana/secondary_metabolism_genes_pathways_with_groups.csv")


    primary_metabolism_v <- distinct(primary_metabolism_0, gene_id) %>%
      merge.data.frame(., all_res, by = "gene_id") %>%
      arrange(RNA_padj)

    secondary_metabolism_v <- distinct(secondary_metabolism_0, gene_id) %>%
      merge.data.frame(., all_res, by = "gene_id") %>%
      arrange(RNA_padj)

    ################
    # Cohen SSE
    sse <- read.table("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/papers/Cohen_14 SSE/gene_list.txt",
      sep = "\t", header = T
    )
    names(sse)[1] <- "gene_id"
    sse$gene_id <- toupper(sse$gene_id)
    sse_merged <- merge.data.frame(sse, all_res, by = "gene_id") %>% dplyr::select(-Title)
    ################
    cat(".")
    # GO term offsprings
    # child_terms_epigenetic <- offspring_fun("GO:0040029") # epigenetic regulation of_gene expression
    child_terms_chromatin_org <- offspring_fun("GO:0006325") # chromatin organization
    child_terms_chromatin_rem <- offspring_fun("GO:0006338") # chromatin remodeling
    cat(".")
    child_terms_defence <- offspring_fun("GO:0006952") # defense response
    child_terms_stress <- offspring_fun("GO:0006950") # response to stress
    cat(".")
    child_terms_biotic <- offspring_fun("GO:0009607") # response to biotic stimulus
    child_terms_abiotic <- offspring_fun("GO:0009628") # response to abiotic stimulus

    ################
    cat(".")
    xl_0_list <- list(
      DNA_methyltransferase = rbind(
        all_res[grep("^2\\.1\\.1\\.37", all_res$EC), ], # EC_2.1.1.37
        all_res[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(all_res$Protein.names)), ], # DNA_C5_MT
        all_res[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(all_res$short_description)), ], # DNA_C5_MT
        all_res[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(all_res$Computational_description)), ],
        all_res[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(all_res$GO.biological.process)), ] # DNA_C5_MT
      ),
      Histone_Lysine_MTs = rbind(
        HLM_mto1, # %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep("SET domain", gsub("SET-domain", "SET domain", all_res$short_description)), ], # %>% dplyr::filter(RNA_pvalue < 0.05), # SET_domain
        # all_res[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(all_res$Gene_description)),], # %>% dplyr::filter(RNA_pvalue < 0.05), # SDG
        all_res[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(all_res$Symbol)), ], # %>% dplyr::filter(RNA_pvalue < 0.05), # SDG
        all_res[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(all_res$old_symbols)), ], # %>% dplyr::filter(RNA_pvalue < 0.05), # SDG
        all_res[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(all_res$Protein.names)), ], # %>% dplyr::filter(RNA_pvalue < 0.05), # SDG
        all_res[grep("class v-like sam-binding methyltransferase", tolower(all_res$Protein.families)), ] # %>% dplyr::filter(RNA_pvalue < 0.05) # SAM_MT
      ),
      RdDM_pathway = rddm_mto1, # %>% dplyr::filter(RNA_pvalue < 0.05),

      Royal_Family_Proteins = RF_mto1, # %>% dplyr::filter(RNA_pvalue < 0.05),
      DNA_deMTs = DNA_deMTs,
      histone_deMTs = histone_deMTs,

      chromatin_remodeling = rbind(
        all_res[grep("chromatin remodeling|chromatin remodeler", tolower(all_res$GO.biological.process)), ], # %>% dplyr::filter(RNA_pvalue < 0.05), # chromatin_remodeling_BP
        all_res[grep("chromatin remodeling|chromatin remodeler", tolower(all_res$Gene_description)), ], # %>% dplyr::filter(RNA_pvalue < 0.05), # chromatin_remodeling_pro_name
        all_res[grep("chromatin remodeling|chromatin remodeler", tolower(all_res$Computational_description)), ], # %>% dplyr::filter(RNA_pvalue < 0.05), # chromatin_remodeling_pro_name
        all_res[grep("chromatin remodeling|chromatin remodeler", tolower(all_res$Function)), ], # %>% dplyr::filter(RNA_pvalue < 0.05), # chromatin_remodeling_pro_name
        all_res[grep("chromatin remodeling|chromatin remodeler", tolower(all_res$short_description)), ] # %>% dplyr::filter(RNA_pvalue < 0.05) # chromatin remodeling_function
      ),
      other_methylation = rbind(
        all_res[grep("class i-like sam-binding methyltransferase", tolower(all_res$Protein.families)), ], # %>% dplyr::filter(RNA_pvalue < 0.05), # SAM_MT
        all_res[grep("class iv-like sam-binding methyltransferase", tolower(all_res$Protein.families)), ] # %>% dplyr::filter(RNA_pvalue < 0.05)#, # SAM_MT
        # all_res[grep("^2\\.1\\.1\\.", all_res$EC.number),] %>% dplyr::filter(RNA_pvalue < 0.05) # EC_2.1.1
      ),
      chromatin_organization_related = rbind(
        all_res[grep("chromosome", tolower(all_res$GO.cellular.component)), ] %>% dplyr::filter(RNA_pvalue < 0.05), # chromosome_CC
        all_res[grep("chromatin", tolower(all_res$GO.cellular.component)), ] %>% dplyr::filter(RNA_pvalue < 0.05), # chromatin_CC
        all_res[grep("chromatin", tolower(all_res$Function)), ] %>% dplyr::filter(RNA_pvalue < 0.05), # chromatin_function
        # all_res[grep("nucleosome", tolower(all_res$GO.cellular.component)), ] %>% dplyr::filter(RNA_pvalue < 0.05), # nucleosome_CC
        # all_res[grep("zinc finger", tolower(all_res$Protein.names)),],  %>% dplyr::filter(RNA_pvalue < 0.05), # zinc_finger_pro_name
        all_res[grep("jumonji", tolower(all_res$short_description)), ] %>% dplyr::filter(RNA_pvalue < 0.05), # jumonji
        all_res[grep("helicase domain", tolower(all_res$short_description)), ] %>% dplyr::filter(RNA_pvalue < 0.05), # helicase_domain
        all_res[grep("histone h3 acetylation", tolower(all_res$GO.biological.process)), ] %>% dplyr::filter(RNA_pvalue < 0.05), # histone_H3_acetylation_BP
        all_res[grep_position(child_terms_chromatin_org), ] %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep_position(child_terms_chromatin_rem), ] %>% dplyr::filter(RNA_pvalue < 0.05)
      ),
      methionine_biosynthesis = all_res[grep("ath00270", all_res$KEGG_pathway), ],
      alan_asp_glut_metabolism = all_res[grep("ath00250", all_res$KEGG_pathway), ],
      gly_ser_and_threo_metabolism = all_res[grep("ath00260", all_res$KEGG_pathway), ],
      lysine_biosynthesis = all_res[grep("ath00300", all_res$KEGG_pathway), ],
      sulfur_biosynthesis = all_res[grep("ath00920", all_res$KEGG_pathway), ],
      primary_metabolism = primary_metabolism_v,
      secondary_metabolism = secondary_metabolism_v,
      AA_transporters = rbind(
        all_res[grep("amino acid transport", tolower(all_res$Short_description)), ] %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep("amino acid transport", tolower(all_res$Computational_description)), ] %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep("amino acid transport", tolower(all_res$Function)), ] %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep("amino acid transport", tolower(all_res$note)), ] %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep("amino acid transport", tolower(all_res$GO.biological.process)), ] %>% dplyr::filter(RNA_pvalue < 0.05)
      ),
      transporters = rbind(
        all_res[grep("transporter", tolower(all_res$Short_description)), ] %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep("transporter", tolower(all_res$Computational_description)), ] %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep("transporter", tolower(all_res$Function)), ] %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep("transporter", tolower(all_res$note)), ] %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep("transport", tolower(all_res$GO.biological.process)), ] %>% dplyr::filter(RNA_pvalue < 0.05),
        all_res[grep("transport", tolower(all_res$GO.molecular.function)), ] %>% dplyr::filter(RNA_pvalue < 0.05)
      ),
      Cohen_SSE = sse_merged, # %>% dplyr::filter(RNA_pvalue < 0.05),

      Ash_SSE = Ash_merged,


      # keep al rows the contain at least one DMR in list of genes (from GO ID)
      # epigenetic_reg._gene_expression = all_res[grep_position(child_terms_epigenetic), ] %>%
      #  dplyr::filter(!(is.na(CG_Genes) & is.na(CHG_Genes) & is.na(CHH_Genes) & is.na(CG_Promoters) & is.na(CHG_Promoters) & is.na(CHH_Promoters))),

      defense_response = all_res[grep_position(child_terms_defence), ] %>%
        dplyr::filter(!(is.na(CG_Genes) & is.na(CHG_Genes) & is.na(CHH_Genes) & is.na(CG_Promoters) & is.na(CHG_Promoters) & is.na(CHH_Promoters))),
      response_to_stress = all_res[grep_position(child_terms_stress), ] %>% dplyr::filter(RNA_pvalue < 0.05),
      # dplyr::filter(!(is.na(CG_Genes) & is.na(CHG_Genes) & is.na(CHH_Genes) & is.na(CG_Promoters) & is.na(CHG_Promoters) & is.na(CHH_Promoters))),

      response_to_biotic = all_res[grep_position(child_terms_biotic), ] %>%
        dplyr::filter(!(is.na(CG_Genes) & is.na(CHG_Genes) & is.na(CHH_Genes) & is.na(CG_Promoters) & is.na(CHG_Promoters) & is.na(CHH_Promoters))),
      response_to_abiotic = all_res[grep_position(child_terms_abiotic), ] %>%
        dplyr::filter(!(is.na(CG_Genes) & is.na(CHG_Genes) & is.na(CHH_Genes) & is.na(CG_Promoters) & is.na(CHG_Promoters) & is.na(CHH_Promoters)))
    )

    cat(".")

    # empty df with correct headers
    tmp_df <- data.frame(matrix(ncol = length(names(xl_0_list[[1]])), nrow = 0))
    names(tmp_df) <- names(xl_0_list[[1]])

    # create data frame to make it unique by rows
    for (i in 1:length(xl_0_list)) {
      tmp_loop_df <- xl_0_list[[i]]
      tmp_loop_df$tmp1 <- names(xl_0_list)[i]

      ### make it unique or not (if each gene will be in one group or all related groups)
      tmp_loop_df$tmp2 <- paste0(tmp_loop_df$gene_id, tmp_loop_df$CG_Genes, tmp_loop_df$CHG_Genes, tmp_loop_df$CHH_Genes, tmp_loop_df$CG_Promoters, tmp_loop_df$CHG_Promoters, tmp_loop_df$CHH_Promoters)
      if (make_it_unique == F) {
        tmp_loop_df <- tmp_loop_df %>%
          distinct(tmp2, .keep_all = TRUE) %>%
          dplyr::select(-tmp2)
      }

      tmp_df <- rbind(tmp_df, tmp_loop_df)
    }

    cat(".")

    if (make_it_unique) {
      tmp_df <- tmp_df %>%
        distinct(tmp2, .keep_all = TRUE) %>%
        dplyr::select(-tmp2)
    }

    # split again to a list
    xl_list <- split(tmp_df, tmp_df$tmp1)
    xl_list <- xl_list[unique(tmp_df$tmp1)] # order it as suppose to be...

    ################################################################################
    cat(".")
    ################################################################################
    remove_dup_DMR <- function(y) {
      y <- as.character(unique(unlist(strsplit(y, ","))))
      paste(y, collapse = ",")
    }


    ################ eddit xl_list data frames
    xl_list <- lapply(xl_list, function(x) {
      # Group by gene_id and summarize (add DMRs values for the gene, make in onw row for each gene)
      x <- x %>%
        group_by(gene_id) %>%
        summarise(
          across(contains("_Genes") | contains("_Promoters"), ~ remove_dup_DMR(paste(., collapse = ","))), # apply remove_dup_DMR function for DMR columns
          across(!contains("_Genes") & !contains("_Promoters"), first) # for other columns
        ) %>%
        as.data.frame() %>%
        mutate(across(contains("_Genes") | contains("_Promoters"), ~ str_replace_all(.x, "NA", ""))) %>%
        mutate(across(contains("_Genes") | contains("_Promoters"), ~ str_replace_all(.x, ",", ", "))) %>% # change delimiter
        # mutate(gene = replace_na(gene, "")) %>%

        relocate(contains("_Genes") | contains("_Promoters"), .before = Symbol) %>%
        # relocate("Function..CC.", .after = "short_description") %>%
        # relocate("gene_model_type", .before = "DOI.ID") %>%
        select(-tmp1) %>%
        arrange(RNA_padj)

      #######

      return(x)
    })

    ################################################################################

    # make only 'epigenetic' sheet unique with others
    dup_genes <- lapply(xl_list[-length(xl_list)], function(x) {
      x[[1]]
    }) %>%
      unlist() %>%
      as.character() %>%
      unique()

    tmp_arg <- xl_list[["epigenetic_reg._gene_expression"]]
    tmp_arg <- tmp_arg[!tmp_arg$gene_id %in% dup_genes, ]
    xl_list[["epigenetic_reg._gene_expression"]] <- tmp_arg

    ################################################################################

    ################ # add beck original columns of 'sse' and 'royal family proteins' and metabolism
    xl_list[["Cohen_SSE"]] <- merge.data.frame(sse, xl_list[["Cohen_SSE"]], by = "gene_id") %>%
      relocate("Title", .after = "Symbol") %>%
      arrange(RNA_padj) %>%
      arrange(Title)

    xl_list[["Royal_Family_Proteins"]] <- merge.data.frame(RF, xl_list[["Royal_Family_Proteins"]], by = "gene_id") %>%
      relocate(c("Agenet.Tudor.Domain", "Other.Domain"), .after = "Symbol") %>%
      arrange(RNA_padj)

    xl_list[["primary_metabolism"]] <- merge.data.frame(primary_metabolism_0, xl_list[["primary_metabolism"]], by = "gene_id") %>%
      relocate(c("metabolite_group", "pathway_names", "kegg_pathways"), .after = "Symbol") %>%
      arrange(RNA_padj) %>%
      arrange(metabolite_group)

    xl_list[["secondary_metabolism"]] <- merge.data.frame(secondary_metabolism_0, xl_list[["secondary_metabolism"]], by = "gene_id") %>%
      relocate(c("metabolite_group", "pathway_names", "kegg_pathways"), .after = "Symbol") %>%
      arrange(RNA_padj) %>%
      arrange(metabolite_group)

    ################################################################################

    # pie plots
    pie_groups <- function(name, is.padj) {
      if (is.padj) {
        x <- xl_list[[name]] %>% filter(RNA_padj < 0.05)
      } else {
        x <- xl_list[[name]] %>% filter(RNA_pvalue < 0.05)
      }

      x_up <- x %>%
        filter(RNA_log2FC > 0) %>%
        nrow()
      x_down <- x %>%
        filter(RNA_log2FC < 0) %>%
        nrow()
      x_total <- x_up + x_down

      pres_up <- round((x_up / x_total) * 100, 1)
      pres_down <- round((x_down / x_total) * 100, 1)

      pie_data <- data.frame(
        group = c(paste0(pres_up, "%"), paste0(pres_down, "%")),
        value = c(x_up, x_down)
      )

      label_up <- paste0(pie_data$group[1], " (", pie_data$value[1], ")")
      label_down <- paste0(pie_data$group[2], " (", pie_data$value[2], ")")

      try({
        if (x_total == 0) {
          plot.new()
          text(0.5, 0.5, "0 DEGs", col = "gray50")
          title(main = paste0("\n\n\n", gsub("_", " ", name)))
          mtext(paste0(rep("_", 24), collapse = ""), side = 1, col = "gray")
        } else {
          pie(c(Up = x_up, Down = x_down),
            labels = c(label_up, label_down),
            main = paste0("\n\n\n", gsub("_", " ", name)),
            col = c("#d96c6c", "#6c96d9"),
            border = "white",
            radius = 0.6
          )
          mtext(paste0(rep("_", 30), collapse = ""), side = 1, col = "gray")
        }
      })
    }

    svg(file = paste0(output_res, "pie_", treatment, "_up_or_down_padj.svg"), width = 9, height = 9, family = "serif")
    par(mfrow = c(5, 5), mar = c(0, 0, 3, 0), oma = c(0, 0, 0, 0))
    for (i_name in names(xl_list)) {
      pie_groups(i_name, is.padj = TRUE)
    }
    dev.off()

    svg(file = paste0(output_res, "pie_", treatment, "_up_or_down_pValue.svg"), width = 9, height = 9, family = "serif")
    par(mfrow = c(5, 5), mar = c(0, 0, 3, 0), oma = c(0, 0, 0, 0))
    for (i_name in names(xl_list)) {
      pie_groups(i_name, is.padj = FALSE)
    }
    dev.off()

    ################################################################################
    clean_ASCII <- function(x) {
      x <- gsub("\001", " ", x)
      x <- gsub("\002", " ", x)
      x <- gsub("\036", " ", x)

      # x = gsub("[[:punct:]]", " ", x)
      # x = iconv(x, from = 'UTF-8', to = 'ASCII')
      return(x)
    }
    ################
    xl_headers <- names(xl_list[[1]])
    numeric_cols <- grep("RNA_|CG_|CHG_|CHH_", xl_headers)
    p_cols <- grep("RNA_p", xl_headers)
    lfc_cols <- grep("RNA_log2FC|CG_|CHG_|CHH_", xl_headers)
    other_cols <- (grep("CHH_Promoters", xl_headers) + 2):length(xl_headers)
    ################
    # save and edit EXCEL
    wb <- createWorkbook()
    # Define styles
    style_up <- createStyle(fontName = "Times New Roman", bgFill = "#f59d98")
    style_down <- createStyle(fontName = "Times New Roman", bgFill = "#c3ccf7")
    style_p <- createStyle(fontName = "Times New Roman", bgFill = "#f7deb0")
    style_other <- createStyle(fontName = "Times New Roman", bgFill = "#daf7d7")
    cell_n_font_style <- createStyle(fontName = "Times New Roman", border = "TopBottomLeftRight", borderColour = "black")
    header_style <- createStyle(fontName = "Times New Roman", textDecoration = "bold", border = "TopBottomLeftRight", borderStyle = "double")

    for (sheet_name in names(xl_list)) {
      addWorksheet(wb, sheet_name)
      #  setColWidths(wb, sheet_name, cols = other_cols, widths = 10)
      #  setColWidths(wb, sheet_name, cols = lfc_cols, widths = 4)
      #  setColWidths(wb, sheet_name, cols = p_cols, widths = 6)

      df <- xl_list[[sheet_name]]
      df <- data.frame(lapply(df, clean_ASCII))
      df[, numeric_cols] <- sapply(df[, numeric_cols], as.numeric)

      writeData(wb, sheet_name, df)

      addStyle(wb, sheet_name, style = cell_n_font_style, rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
      addStyle(wb, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)

      conditionalFormatting(wb, sheet_name, cols = p_cols[1], rows = 2:(nrow(df) + 1), rule = "<0.05", style = style_p)
      conditionalFormatting(wb, sheet_name, cols = p_cols[2], rows = 2:(nrow(df) + 1), rule = "<0.05", style = style_p)

      for (col in lfc_cols) {
        conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df) + 1), rule = ">0", style = style_up)
        conditionalFormatting(wb, sheet_name, cols = col, rows = 2:(nrow(df) + 1), rule = "<0", style = style_down)
      }

      for (col in other_cols[other_cols %% 2 == 0]) {
        conditionalFormatting(wb, sheet_name, style = style_other, rule = "!=0", rows = 2:(nrow(df) + 1), cols = col, gridExpand = TRUE)
        conditionalFormatting(wb, sheet_name, style = style_other, rule = "==0", rows = 2:(nrow(df) + 1), cols = col, gridExpand = TRUE)
      }
      cat(".")
    }
    saveWorkbook(wb, paste0(output_res, treatment, unique_or_not, "_groups.xlsx"), overwrite = T)
  }
  cat(" done\n")
}
