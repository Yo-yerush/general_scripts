###########################
#                         #
#  $ conda activate Renv  #
#                         #
###########################

output_dir <- "/home/yoyerush/yo/methylome_pipeline/mto3_delta_df/"

methyl_files_type <- "CX_report"

wt_name <- "wt"
mutant_name <- "cmt3"

wt_path <- c(
    "/home/yoyerush/yo/methylome_pipeline/other_mutants/stroud_et_al_2013/bismark_results/wt_2_bismark_se.CX_report.txt.gz",
    "/home/yoyerush/yo/methylome_pipeline/other_mutants/stroud_et_al_2013/bismark_results/wt_3_bismark_se.CX_report.txt.gz"
)
mutant_path <- c(
    "/home/yoyerush/yo/methylome_pipeline/other_mutants/stroud_et_al_2013/bismark_results/cmt3_bismark_se.CX_report.txt.gz"
)

n.cores <- 6

############# ############# ############# #############

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(DMRcaller)
library(parallel)
library(data.table)

source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/trimm_and_rename_seq.R")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.At/main/scripts/load_replicates.R")

dir.create(output_dir, showWarnings = F)
############# ############# ############# #############

# load replicates data
var_args <- list(
    list(path = wt_path, name = wt_name),
    list(path = mutant_path, name = mutant_name)
)

n.cores.load <- ifelse(n.cores > 1, 2, 1)
load_vars <- mclapply(var_args, function(x) load_replicates(x$path, n.cores, x$name, T, methyl_files_type), mc.cores = n.cores.load)

meth_wt <- trimm_and_rename(load_vars[[1]]) %>%
    as.data.frame() %>%
    mutate(Proportion = readsM / readsN)

meth_mutant <- trimm_and_rename(load_vars[[2]]) %>%
    as.data.frame() %>%
    mutate(Proportion = readsM / readsN)

############# ############# ############# #############

meth_delta <- meth_wt %>%
    mutate(delta = meth_mutant$Proportion - meth_wt$Proportion) %>%
    select(seqnames, start, end, delta, context) %>%
    filter(!is.nan(delta))


CG_delta <- meth_delta[meth_delta$context == "CG", ]
CHG_delta <- meth_delta[meth_delta$context == "CHG", ]
CHH_delta <- meth_delta[meth_delta$context == "CHH", ]

write.csv(rbind(CG_delta, CHG_delta, CHH_delta),
    paste0(output_dir, mutant_name, "_delta_df.csv.gz"),
    row.names = FALSE,
    quote = FALSE
)
#
