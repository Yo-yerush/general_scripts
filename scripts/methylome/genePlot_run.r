library(DMRcaller)
library(GenomicFeatures)
library(dplyr)
library(parallel)

tair_id <- "AT3G01120"

#methylome_at_results <- NULL
methylome_at_results <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_240325/Methylome.At/results/mto1_vs_wt"

output_path <- "/home/yoyerush/yo/genePlot_080725"

# upload CX files
var1_path <- c(
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S18/methylation_extractor/S18_R1_bismark_bt2_pe.CX_report.txt",
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S19/methylation_extractor/S19_R1_bismark_bt2_pe.CX_report.txt"
)

var2_path <- c(
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S9/methylation_extractor/S9_R1_bismark_bt2_pe.CX_report.txt",
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S10/methylation_extractor/S10_R1_bismark_bt2_pe.CX_report.txt",
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S11/methylation_extractor/S11_R1_bismark_bt2_pe.CX_report.txt"
)



###################################

var1_name <- "WT"
var2_name <- "mto1"

if (n_cores >= length(c(var1_path, var2_path))) {
  ncores_par <- n_cores
} else {
  ncores_par <- 1
}

vars_files <- mclapply(c(var1_path, var2_path), readBismark, mc.cores = ncores_par)
var1 <- vars_files[1:length(var1_path)]
var2 <- vars_files[1:length(var2_path)]

var1_pool <- poolMethylationDatasets(GRangesList(var1))
var2_pool <- poolMethylationDatasets(GRangesList(var2))

var1_pool <- rename_seq(var1_pool)
var2_pool <- rename_seq(var2_pool)

###################################

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/methylome/genePlot_script.r")
genePlot_fun(tair_id, var1_pool, var2_pool, var1_name, var2_name, methylome_at_results, output_path)
