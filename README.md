-----------------------------------------------------------------
-----------------------------------------------------------------

# WGBS pipeline

#### * Download fastq files from SRA and run [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/#:~:text=Bismark%20is%20a%20program%20to%20map%20bisulfite%20treated,the%20methylation%20levels%20of%20their%20samples%20straight%20away.) for *dml3* samples ([Zhejiang University](https://www.ncbi.nlm.nih.gov/sra/SRX4698864))
```
./download_fastq_from_SRA_short.sh "SRR7848067 SRR7848068 SRR7848069 SRR7848070" "/PATH/TO/dml3_NS"
```

#### * Create a sample table file (*tab* delimiter)
```
dml3_1    PATH/TO/FILE/dml3_1_R1.fastq
dml3_1    PATH/TO/FILE/dml3_1_R2.fastq
dml3_2    PATH/TO/FILE/dml3_2_R1.fastq
dml3_2    PATH/TO/FILE/dml3_2_R2.fastq
wt_1    PATH/TO/FILE/wt1_R1.fastq
wt_1    PATH/TO/FILE/wt1_R2.fastq
wt_2    PATH/TO/FILE/wt2_R1.fastq
wt_2    PATH/TO/FILE/wt2_R2.fastq
```
#### * Run Bismark to get only '**.CX_report.txt**' file
can run without '--cx' option to get all output files
```
.run_bismark_yo.sh -s ./dml3_NS/bismark_dml3_samples.txt -g /PATH/TO/TAIR10_chr_all.fas.gz -o ./dml3_NS/bismark_results -n 32 -m 16G --cx
```

#### * Use [Methylome.At](https://github.com/Yo-yerush/Methylome.At) downstream pipeline
```

```

-----------------------------------------------------------------
-----------------------------------------------------------------

# Calculate and plot '*delta*' methylation levels
*[Stroud et al. 2013](https://pubmed.ncbi.nlm.nih.gov/23313553/). (SRA experiment: [SRP014726](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP014726&o=biosample_s%3Aa%3Bacc_s%3Aa))*
#### Download '[*Stroud et al. (2013)*](https://pubmed.ncbi.nlm.nih.gov/23313553/)' SRA files  (*SRA experiment: [SRP014726](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP014726&o=biosample_s%3Aa%3Bacc_s%3Aa)*) and use [mutants compare_delta_df.r](https://github.com/Yo-yerush/general_scripts/blob/main/delta_df_from_wig_script.r) script
This will save '**.csv**' files of the total-methylation delta (mutants compared to WT)
in this example, *mto1* mutant '.csv' file created by '**.CX_report.txt**' file, using [mutants compare_delta_df.r](https://github.com/Yo-yerush/general_scripts/blob/main/delta_df_from_CX_report_script.r) script.

#### * Then to create a ChrPlot use the following script
```
library(ggplot2)
library(dplyr)
library(GenomicRanges)

output_dir <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/mutants_figs"

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/windowSize_for_GRanges_mcol.r")
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/ChrPlots_CX_yo.R")

read_mut_file <- function(mut_name) {
    cat(paste0("read ", mut_name, " file..."))
    x <- read.csv(paste0("PATH/TO/mutants_figs/", mut_name, "_delta_df.csv.gz"))

    return(list(
        cg = x[x$context == "CG", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta"),
        chg = x[x$context == "CHG", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta"),
        chh = x[x$context == "CHH", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta")
    ))
    cat(" done\n")
}

mto1 <- read_mut_file("mto1")
met1 <- read_mut_file("met1")
cmt2 <- read_mut_file("cmt2")
cmt3 <- read_mut_file("cmt3")
ddm1 <- read_mut_file("ddm1")


svg(paste0(output_dir, "/ChrPlot_test_stroud_all.svg"), width = 7, height = 4, family = "serif")
ChrPlots_CX_all(
    meth_var_list = "test_stroud_all",
    meth_names = list(mto1,met1,cmt2,cmt3,ddm1), c("mto1","met1","cmt2","cmt3","ddm1"),
    y_max_cg = 0,
    y_max_chg = 0.2,
    y_max_chh = 0.05,
    y_mid_cg = NULL,
    y_mid_chg = 0,
    y_mid_chh = 0,
    y_min_cg = -1,
    y_min_chg = -0.5,
    y_min_chh = -0.15,
    italic_legend_names = TRUE,
    ylab_suffix = "(Δ)",
    y_title_cex = 1.1,
    TE_as_gr = "tair10"
    )
dev.off()
```
#### ChrPlots:
![fig](https://github.com/Yo-yerush/general_scripts/blob/main/ChrPlot_test_stroud_all.svg)

-----------------------------------------------------------------
-----------------------------------------------------------------

