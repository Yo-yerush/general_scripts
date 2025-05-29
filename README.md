### Download fastq files from SRA and run Bismark for dml3 samples
```
./download_fastq_from_SRA_short.sh "SRR7848067 SRR7848068 SRR7848069 SRR7848070" "/home/yoyerush/yo/methylome_pipeline/other_mutants/dml3_NS"
```

### Create a sample table file (example):
```
dml_1    PATH/TO/FILE/mt1_R1.fastq
dml_1    PATH/TO/FILE/mt1_R2.fastq
dml_2    PATH/TO/FILE/mt2_R1.fastq
dml_2    PATH/TO/FILE/mt2_R2.fastq
wt_1    PATH/TO/FILE/wt1_R1.fastq
wt_1    PATH/TO/FILE/wt1_R2.fastq
wt_2    PATH/TO/FILE/wt2_R1.fastq
wt_2    PATH/TO/FILE/wt2_R2.fastq
```
### run BISMARK to get just 'CX_report' file (can run without --cx option to get all ouput files)
```
.run_bismark_yo.sh -s dml3_NS/bismark_dml3_samples.txt -g /home/yoyerush/TAIR10_chr_all.fas.gz -o ./dml3_NS/bismark_results -n 32 -m 16G --cx
```

##################################################
# Calculate 'delta' methylation levels from '.wig' files ('stroud et al. 2013' example)
### Run [mutants compare_delta_df.r](https://github.com/Yo-yerush/general_scripts/blob/main/mutants_compare_delta_df.r) script

### Than to create a ChrPlot:
```
library(ggplot2)
library(dplyr)
library(GenomicRanges)

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/windowSize_for_GRanges_mcol.r")
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/ChrPlots_CX_from_GRanges.R")

read_mut_file <- function(mut_name) {
    x <- read.csv(paste0("PATH/TO/mutants_figs/", mut_name, "_delta_df.csv.gz"))

    return(list(
        cg = x[x$context == "CG", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta"),
        chg = x[x$context == "CHG", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta"),
        chh = x[x$context == "CHH", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta")
    ))
}

met1 <- read_mut_file("met1")
cmt2 <- read_mut_file("cmt2")
cmt3 <- read_mut_file("cmt3")

ChrPlots_CX("MTs_stroud", list(met1$cg,cmt2$cg,cmt3$cg), c("met1","cmt2","cmt3"), "CG", y_max = 0.1, y_mid = 0, y_min = -1, output_dir="PATH/TO/mutants_figs")
ChrPlots_CX("MTs_stroud", list(met1$cg,cmt2$cg,cmt3$cg), c("met1","cmt2","cmt3"), "CHG", y_max = 0.2, y_mid = 0, y_min = -0.5, output_dir="PATH/TO/mutants_figs")
ChrPlots_CX("MTs_stroud", list(met1$cg,cmt2$cg,cmt3$cg), c("met1","cmt2","cmt3"), "CHH", y_max = 0.05, y_mid = 0, y_min = -0.15, output_dir="PATH/TO/mutants_figs")
```
