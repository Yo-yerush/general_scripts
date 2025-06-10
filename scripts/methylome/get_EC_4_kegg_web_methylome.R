library(dplyr)
library(tidyr)
library(org.At.tair.db)
library(KEGGREST)
library(pathview)
library(purrr)

# by_TPM = T # if TRUE the pValue is raw TPM calculation. if FALSE the pValue is DEseq2 output values

map_ids <- c(
  "03020",
  "04122",
  "04712",
  "00270",
  "00300",
  "00310",
  "00340",
  "00920",
  "00020",
  "00500",
  "00520",
  "00260",
  "00290",
  "00330",
  "00250",
  "00620" # ,
  # "03082", # cant find tairs by KEGGREST and org.At.tair.db
  # "03083", # cant find tairs by KEGGREST and org.At.tair.db
  # "01230", # cant find tairs by KEGGREST and org.At.tair.db
  # "01232", # cant find tairs by KEGGREST and org.At.tair.db
)


for (treatment in c("mto1", "mto3", "dCGS", "SSE_high", "SSE_low", "SSE_high_vs_SSE_low")) {
  for (context in c("CG", "CHG", "CHH", "all")) {
    for (annotation in c("Genes", "Promoters")) {
      for (n.map in map_ids) {
        for (is.EC in c(TRUE)) { # FALSE
          try({
            if (treatment == "SSE_high_vs_SSE_low") {
              treatment_name <- treatment
            } else {
              treatment_name <- ifelse(grepl("mto", treatment), paste0(treatment, "_vs_wt"), paste0(treatment, "_vs_EV"))
            }

            # pathway name for output file
            pathway_info <- keggGet(paste0("map", n.map))
            pathway_name <- paste0(pathway_info[[1]]$NAME, "_ath", n.map)

            # all tairs for ath map
            tairs <- data.frame(gene_id = as.list(org.At.tairPATH2TAIR)[[n.map]])

            # RNA results
            if (context == "all") {
              meth_file <- rbind(
                read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/", treatment_name, "/genome_annotation/CG/", annotation, "_CG_genom_annotations.csv"))[, c("gene_id", "Symbol", "log2FC", "pValue")],
                read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/", treatment_name, "/genome_annotation/CHG/", annotation, "_CHG_genom_annotations.csv"))[, c("gene_id", "Symbol", "log2FC", "pValue")],
                read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/", treatment_name, "/genome_annotation/CHH/", annotation, "_CHH_genom_annotations.csv"))[, c("gene_id", "Symbol", "log2FC", "pValue")]
              )
            } else {
              meth_file <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/", treatment_name, "/genome_annotation/", context, "/", annotation, "_", context, "_genom_annotations.csv"))[, c("gene_id", "Symbol", "log2FC", "pValue")]
            }

            # all EC for the same map
            if (is.EC) {
              ec_ids <- keggLink("enzyme", paste0("map", n.map))
              is.EC <- ifelse(length(ec_ids) != 0, T, F)
            }

            # name_if_ec <- ifelse(is.EC, "_with_EC", "")
            # dir_if_ec <- ifelse(is.EC, "with_EC/", "")
            name_if_ec <- "with_EC_"
            dir_if_ec <- "with_EC/"

            if (is.EC) {
              ec2tair_0 <- keggLink("ath", ec_ids)
              ec2tair <- data.frame(gene_id = ec2tair_0, ec = names(ec2tair_0))
              ec2tair$gene_id <- gsub("ath:", "", ec2tair$gene_id)
              ec2tair$ec <- gsub("ec:", "", ec2tair$ec)

              # merge TAIR and EC ids
              ec_map <- merge.data.frame(tairs, ec2tair, by = "gene_id", all.y = T)

              # group by EC for KEGG ath map (put the fitst tair for web map location)
              ec_map_grouped <- ec_map %>%
                group_by(ec) %>%
                mutate(group_id = dplyr::first(gene_id)) %>%
                ungroup() %>%
                as.data.frame()
              ec_map_grouped <- ec_map_grouped[!duplicated(ec_map_grouped$ec), 2:3]


              merged_tair <- merge.data.frame(ec_map, meth_file, by = "gene_id") %>%
                dplyr::relocate(Symbol, .after = gene_id) %>%
                arrange(pValue) # %>% na.omit()
            } else {
              merged_tair <- merge.data.frame(tairs, meth_file, by = "gene_id") %>%
                arrange(pValue) # %>% na.omit()
              # as.list(org.At.tairENZYME2TAIR)
            }

            # create directories
            path_name <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/KEGG_pathway/"
            dir.create(path_name, showWarnings = F)
            dir.create(paste0(path_name, dir_if_ec), showWarnings = F)
            dir.create(paste0(path_name, dir_if_ec, treatment_name), showWarnings = F)
            dir.create(paste0(path_name, dir_if_ec, treatment_name, "/", annotation), showWarnings = F)
            path_save <- paste0(path_name, dir_if_ec, treatment_name, "/", annotation, "/", context, "/")
            dir.create(path_save, showWarnings = F)
            dir.create(paste0(path_save, "pathview"), showWarnings = F)

            # save RNA seq results for KEGG map
            write.csv(merged_tair, paste0(path_save, treatment, "_values_", pathway_name, "_", context, "_", annotation, name_if_ec, ".csv"), row.names = F)

            # save pathway view map (KEGG)
            if (is.EC) {
              # group by EC
              pathview_df <- merged_tair %>%
                filter(pValue < 0.05) %>%
                distinct(gene_id, .keep_all = T) %>%
                group_by(ec) %>%
                summarise(
                  gene_id = unique(gene_id)[1],
                  # Symbol  = paste(unique(Symbol), collapse = ";"),
                  log2FCs = list(log2FC) # collect each group's log2FC into a list
                ) %>%
                tidyr::unnest_wider(log2FCs, names_sep = "_") %>%
                # fillo into non-NA values for multiple colors
                rowwise() %>%
                mutate(
                  filled = list({
                    x <- c_across(starts_with("log2FCs_"))
                    vals <- x[!is.na(x)]
                    n <- length(x)
                    rep(vals, each = n / length(vals))
                  })
                ) %>%
                ungroup() %>%
                dplyr::select(-starts_with("log2FCs_")) %>%
                tidyr::unnest_wider(filled, names_sep = "") %>%
                as.data.frame() %>%
                dplyr::select(-ec)
            } else {
              pathview_df <- merged_tair %>%
                filter(pValue < 0.05) %>%
                distinct(gene_id, .keep_all = T) %>%
                dplyr::select(-Symbol, -pValue)
            }

            row.names(pathview_df) <- pathview_df$gene_id
            pathview_df <- dplyr::select(pathview_df, , -gene_id)

            setwd(paste0(path_save, "pathview"))
            pathview(
              gene.data = pathview_df,
              # cpd.data = sim.cpd.data,
              pathway.id = n.map,
              gene.idtype = "kegg",
              species = "ath",
              out.suffix = paste0(treatment, "_", pathway_name, "_", context, "_", annotation, name_if_ec),
              keys.align = "y",
              kegg.native = T,
              key.pos = "topright",
              sign.pos = "bottomleft",
              cpd.lab.offset = -1,
              low = list(gene = "#4949f5", cpd = "#47e047"),
              mid = list(gene = "gray95", cpd = "gray95"),
              high = list(gene = "#f14e4e", cpd = "#f0d851"),
              dsicrete = T,
              # limit = c(min(pathview_df$log2FC), max(pathview_df$log2FC)),
              limit = 2,
              bins = 20,
              match.data = T,
              multi.state = T,
              same.layer = F
            )
            # Remove the file named "yo" if it exists
            file.remove(paste0("ath", n.map, ".xml"))
            file.remove(paste0("ath", n.map, ".png"))
            setwd("~")

            # transform Inf to max values
            merged_tair$log2FC[merged_tair$log2FC == Inf] <- unique(sort(merged_tair$log2FC, decreasing = T))[2]
            merged_tair$log2FC[merged_tair$log2FC == -Inf] <- unique(sort(merged_tair$log2FC, decreasing = F))[2]

            if (is.EC) {
              rmv_pos <- grep("gene_id|pValue", names(merged_tair))
            } else {
              rmv_pos <- grep("pValue", names(merged_tair))
            }


            # merged_up = merged_tair[merged_tair$pValue < 0.05 & merged_tair$log2FC > 0, -rmv_pos] %>% arrange(desc(log2FC))
            # merged_down = merged_tair[merged_tair$pValue < 0.05 & merged_tair$log2FC < 0, -rmv_pos] %>% arrange(desc(log2FC))
            merged_sig <- merged_tair[merged_tair$pValue < 0.05, -rmv_pos] %>% arrange(desc(log2FC))
            merged_non_sig <- merged_tair[merged_tair$pValue > 0.05, -rmv_pos]
            if (nrow(merged_non_sig) != 0) {
              merged_non_sig$log2FC <- "#fae7b4"
            }

            # transform max max value to 1 and min value to -1
            up_term <- merged_sig$log2FC > 0
            down_term <- merged_sig$log2FC < 0
            merged_sig$log2FC[up_term] <- merged_sig$log2FC[up_term] / max(merged_sig$log2FC[up_term])
            merged_sig$log2FC[down_term] <- -(merged_sig$log2FC[down_term] / min(merged_sig$log2FC[down_term]))

            merged_sig$log2FC[up_term & merged_sig$log2FC < 0.35] <- 0.35
            merged_sig$log2FC[down_term & merged_sig$log2FC > -0.35] <- -0.35

            # group by EC to have  a column of log2FC results (with " " delimiter)
            df_2_save <- rbind(merged_sig, merged_non_sig) %>%
              group_by(.data[[names(merged_sig)[1]]]) %>%
              summarise(log2FC = paste(log2FC, collapse = " ")) %>%
              ungroup() %>%
              as.data.frame()

            if (is.EC) {
              # merge with ec forups for tair ids column (ath map need tair id and not ec)
              df_2_save <- merge.data.frame(ec_map_grouped, df_2_save, by = "ec")[, -1]
            }

            df_2_save$log2FC <- gsub(" #fae7b4", "", df_2_save$log2FC)

            write.table(df_2_save, paste0(path_save, treatment, "_", pathway_name, "_", context, name_if_ec, ".txt"), row.names = F, col.names = F, sep = " ", quote = F)
          })
        }
      }
    }
  }
}
