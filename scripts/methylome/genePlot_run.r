#############3333
# test



# variable names
var1_name <- "WT" # control
var2_name <- "mto1"

tair_id <- "AT3G01120"

# methylome_at_results <- NULL
methylome_at_results <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At_240325/Methylome.At/results/mto1_vs_wt"

output_path <- "/home/yoyerush/yo/genePlot_080725"

# CX files
var1_path <- c(
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S18/methylation_extractor/S18_R1_bismark_bt2_pe.CX_report.txt",
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S19/methylation_extractor/S19_R1_bismark_bt2_pe.CX_report.txt"
)

var2_path <- c(
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S9/methylation_extractor/S9_R1_bismark_bt2_pe.CX_report.txt",
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S10/methylation_extractor/S10_R1_bismark_bt2_pe.CX_report.txt",
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S11/methylation_extractor/S11_R1_bismark_bt2_pe.CX_report.txt"
)

n_cores <- length(c(var1_path, var2_path))

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/methylome/genePlot_script.r")
genePlot_fun(tair_id, var1_path, var2_path, var1_name, var2_name, methylome_at_results, output_path, n_cores, create_legend = TRUE)
