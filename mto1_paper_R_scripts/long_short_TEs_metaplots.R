# run it on darwin (linux server)
library(dplyr)
library(ggplot2)
library(DMRcaller)
library(org.At.tair.db)
library(GenomicFeatures)
library(plyranges)
library(parallel)
library(data.table)

source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/Genes_metaPlot_fun.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/trimm_and_rename_seq.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/load_replicates.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/edit_TE_file.R")

TE_file_path <- "https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/annotation_files/TAIR10_Transposable_Elements.txt"
samples_path_df <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/Methylome.At/samples_table/samples_table_mto1.txt"
n.cores <- 30

#############################################

# load replicates raw 'CX' data
var_table <- read.csv(samples_path_df, header = F, sep = "\t")
vars_vector <- unique(var_table[, 1])
var1_path <- var_table[grep(vars_vector[1], var_table[, 1]), 2]
var2_path <- var_table[grep(vars_vector[2], var_table[, 1]), 2]
var_args <- list(
    list(path = var1_path, name = "wt"),
    list(path = var2_path, name = "mto1")
)
load_vars <- mclapply(var_args, function(x) load_replicates(x$path, n.cores, x$name, T, "CX_report"), mc.cores = 5)
meth_wt <- trimm_and_rename(load_vars[[1]])
meth_mto1 <- trimm_and_rename(load_vars[[2]])

#############################################

TE_df <- read.csv(TE_file_path, sep = "\t")
te_width <- data.frame(
    te_id = TE_df$Transposon_Name,
    width = (TE_df$Transposon_max_End - TE_df$Transposon_min_Start)
)
long_tes <- te_width[te_width$width >= 4000, 1]
short_tes <- te_width[te_width$width <= 500, 1] # %>% sample(10000) # random IDs

#############################################

metaPlot_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots"
long_te_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots/long_TEs"
short_te_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots/short_TEs"
dir.create(metaPlot_path, showWarnings = F)
dir.create(long_te_path, showWarnings = F)
dir.create(short_te_path, showWarnings = F)

setwd(long_te_path)
Genes_metaPlot(meth_wt, meth_mto1, "wt", "mto1", edit_TE_file(TE_df), long_tes, 6, n.cores, is_TE = T)

setwd(short_te_path)
Genes_metaPlot(meth_wt, meth_mto1, "wt", "mto1", edit_TE_file(TE_df), short_tes, 6, n.cores, is_TE = T)

########################
## all TEs run
all_tes <- te_width[, 1]
all_te_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/mto1_long_short_TE_metaplots/all_TEs"
dir.create(all_te_path, showWarnings = F)

setwd(all_te_path)
Genes_metaPlot(meth_wt, meth_mto1, "wt", "mto1", edit_TE_file(TE_df), all_tes, 6, n.cores, is_TE = T)

########################
setwd(metaPlot_path)

# Write summary information to a text file - thanks GPT
summary_info <- paste0(
    "Transposable Elements (TEs) Methylation Analysis Summary\n",
    "======================================================\n\n",
    "Analysis Date: ", Sys.Date(), "\n",
    "Analysis Time: ", Sys.time(), "\n\n",
    "Data Sources:\n",
    "- TE annotation file: ", TE_file_path, "\n",
    "- Samples table: ", samples_path_df, "\n\n",
    "Sample Information:\n",
    "- Wild type (wt): ", length(var1_path), " replicates\n",
    "- mto1 mutant: ", length(var2_path), " replicates\n\n",
    "TE Categories Analyzed:\n",
    "- Total TEs in dataset: ", nrow(TE_df), "\n",
    "- Long TEs (≥4000 bp): ", length(long_tes), " TEs\n",
    "- Short TEs (≤500 bp): ", length(short_tes), " TEs\n",
    "- All TEs analyzed: ", length(all_tes), " TEs\n\n",
    "Analysis Parameters:\n",
    "- Number of cores used: ", n.cores, "\n",
    "- Bin size for metaplots: 6\n",
    "- Context analyzed: CX (all cytosines)\n\n",
    "Output Directories:\n",
    "- Main output: ", metaPlot_path, "\n",
    "- Long TEs: ", long_te_path, "\n",
    "- Short TEs: ", short_te_path, "\n",
    "- All TEs: ", all_te_path, "\n"
)
writeLines(summary_info, file.path(metaPlot_path, "analysis_summary.txt"))