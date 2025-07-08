
tair_id <- "AT3G01120"
expitiment_name <- ""
methylome_at_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At"
methylome_at_results <- "/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/mto1_vs_wt"

output_path <- "/home/yoyerush/yo/genePlot_101124"

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

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/methylome/genePlot_script.r")
genePlot_fun(tair_id, var1_path, var2_path, methylome_at_annotations, methylome_at_results = NULL)