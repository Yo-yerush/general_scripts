Workflows:

1. Download SRA files (as '*.fastq*')
2. WGBS pipeline - using [Bismark software](https://www.bioinformatics.babraham.ac.uk/projects/bismark/#:~:text=Bismark%20is%20a%20program%20to%20map%20bisulfite%20treated,the%20methylation%20levels%20of%20their%20samples%20straight%20away.)
3. [Methylome.At](https://github.com/Yo-yerush/Methylome.At) downstream pipeline

4. Plot **delta** methylation levels from '*.wig*' or '*.CX_report.txt*' files
5. genePlots
6. long/short TEs
7. sub-context analysis

8. RNAseq pipeline - using [RSEM software](https://github.com/deweylab/RSEM)
9. RNAseq downstream pipeline

10. DMRs-DEGs correlations
11. Create 'expression + DMRs' integrated tables from WGBS and RNAseq data
12. Classify into gene groups (epigenetic, metabolism, stress-related, etc.)
13. Methylated TEs near DEGs

14. KEGG pathways - enrichment and viewPlots
15. GO analysis
16. REVIGO scripts

-----------------------------------------------------------------
-----------------------------------------------------------------

## Download fastq files from SRA
#### Download fastq files from SRA and run [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/#:~:text=Bismark%20is%20a%20program%20to%20map%20bisulfite%20treated,the%20methylation%20levels%20of%20their%20samples%20straight%20away.) for *dml3* samples ([Zhejiang University](https://www.ncbi.nlm.nih.gov/sra/SRX4698864))
```bash
./download_fastq_from_SRA_short.sh "SRR7848067 SRR7848068 SRR7848069 SRR7848070" "/PATH/TO/dml3_NS"
```
-----------------------------------------------------------------
-----------------------------------------------------------------

## WGBS pipeline
#### Run [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/#:~:text=Bismark%20is%20a%20program%20to%20map%20bisulfite%20treated,the%20methylation%20levels%20of%20their%20samples%20straight%20away.) for *dml3* samples ([Zhejiang University](https://www.ncbi.nlm.nih.gov/sra/SRX4698864))

 Download *Arabidopsis* reference genome ([TAIR10](https://www.arabidopsis.org/))
 ```bash
 cd /PATH/TO
 wget -O TAIR10_chr_all.fas.gz https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
 ```

Create a sample table file (*tab* delimiter)
```txt
dml3_1    PATH/TO/FILE/dml3_1_R1.fastq
dml3_1    PATH/TO/FILE/dml3_1_R2.fastq
dml3_2    PATH/TO/FILE/dml3_2_R1.fastq
dml3_2    PATH/TO/FILE/dml3_2_R2.fastq
wt_1    PATH/TO/FILE/wt1_R1.fastq
wt_1    PATH/TO/FILE/wt1_R2.fastq
wt_2    PATH/TO/FILE/wt2_R1.fastq
wt_2    PATH/TO/FILE/wt2_R2.fastq
```

Run Bismark to get only '**.CX_report.txt**' file
* *run without '--cx' option to get all output files*
```bash
.run_bismark_yo.sh -s ./dml3_NS/bismark_dml3_samples.txt -g /PATH/TO/TAIR10_chr_all.fas.gz -o ./dml3_NS/bismark_results -n 32 -m 16G --cx
```

#### Use [Methylome.At](https://github.com/Yo-yerush/Methylome.At) downstream pipeline
```
git clone https://github.com/Yo-yerush/Methylome.At.git
cd ./Methylome.At
chmod +x ./setup_env.sh
./setup_env.sh

# create sample table with two columns (tab delimiter):
# SAMPLE PATH

# run pipelines
./Methylome.At.sh PATH/TO/samples_table.txt
./Methylome.At_metaPlots.sh PATH/TO/samples_table.txt
```

-----------------------------------------------------------------
-----------------------------------------------------------------

## Calculate and plot '*delta*' methylation levels
Download '[*Stroud et al. (2013)*](https://pubmed.ncbi.nlm.nih.gov/23313553/)' '**.wig**' files  (*SRA experiment: '[SRP014726](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP014726&o=biosample_s%3Aa%3Bacc_s%3Aa)'*) and use '[mutants compare_delta_df.r](https://github.com/Yo-yerush/general_scripts/blob/main/delta_df_from_wig_script.r)' script.
This will save '**.csv**' files of the total-methylation delta (mutants compared to WT).
In this example, *mto1* mutant '**.csv**' file created by '**.CX_report.txt**' (output from Bismark).

#### Then to create **ChrPlots** run the following commands
*  *use <TE_as_gr = NULL> argument to remove TE density from the plot*
```r
library(ggplot2)
library(dplyr)
library(GenomicRanges)

output_dir <- "PATH/TO/mutants_figs"

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
    meth_var_list = list(mto1,met1,cmt2,cmt3,ddm1),
    meth_names = c("mto1","met1","cmt2","cmt3","ddm1"),
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
#### example output
![fig](https://github.com/Yo-yerush/general_scripts/blob/main/ChrPlot_test_stroud_all.svg)

-----------------------------------------------------------------
-----------------------------------------------------------------

## genePlots of methylations over the gene body and its 2kb upstream region
#### One or multiple TAIR ID(s)
#### This script use 'CX_report' files as input
* *can add DMRs results (using 'Methylome.At' pipeline results)*
#### Run in multiple cores are avilable for CX_repots reading files. in default: n_cores => total samples number 
#### To create **genePlots** use the following script
```r
# DMRcaller package-based functions
# Catoni, Marco, Tsang, MF J, Greco, P A, Zabet, Radu N (2018). “DMRcaller: a versatile R/Bioconductor package for detection and visualization of differentially methylated regions in CpG and non-CpG contexts.” Nucleic Acids Research. doi:10.1093/nar/gky602.
# https://www.bioconductor.org/packages/release/bioc/html/DMRcaller.html
# ---------------------------------------------------------------------#

# TAIR ID(s) to plot
tair_id <- c("AT3G01120", "AT3G17390", "AT1G53480")

# variable names
var1_name <- "WT" # control
var2_name <- "mto1"

# output path (will create if not exist)
output_path <- "/PATH/TO/output_directory"

# CX files path
var1_path <- c(
  "/PATH/TO/wt_1_pe.CX_report.txt",
  "/PATH/TO/wt_1_pe.CX_report.txt"
)

var2_path <- c(
  "/PATH/TO/mto1_1_pe.CX_report.txt",
  "/PATH/TO/mto1_2_pe.CX_report.txt",
  "/PATH/TO/mto1_3_pe.CX_report.txt"
)

n_cores <- length(c(var1_path, var2_path))

# DMRs results
# use 'NULL' if there is no 'Methylome.At' results output
methylome_at_results <- "/PATH/TO/Methylome.At/results/mto1_vs_wt"

# run genePlots scripts
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/methylome/genePlot_script.r")
genePlot_fun(tair_id, var1_path, var2_path, var1_name, var2_name, methylome_at_results, output_path, n_cores, create_legend = TRUE)
```
#### example output

<img
  src="https://github.com/Yo-yerush/general_scripts/blob/main/genePlot_AT3G01120.svg"
  alt="fig"
  width="50%"
/>

-----------------------------------------------------------------
-----------------------------------------------------------------

## DMRs-DEGs correlations
First, make average methylation level for each sample, in each context, and for each gene feature (promoter, CDS, etc.)
* *note: 'average_methylation_over_gene_features' function will run **42** jobs in parallel (max)*
```r
sample_file_path <- "/PATH/TO/samples_table.txt"
output_path <- "/PATH/TO/OUTPUT/DIRECTORY/"

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/methylome/RNAseq_correlation/sample_ave_methylation_levels_from_gene_features.R")
average_methylation_over_gene_features(sample_file_path, output_path)
```

Create a linear plot using **Negative Binomial Regression**.
need to load:
* average methylation over gene features
* DMRs results deirectory (if run Methylome.At pipeline)
* RNA-seq results directoey (both statistic results and norm.Count files). *better use RNAseq downstream pipeline*
```r
var1_name <- "wt"
var2_name <- "SSE_high"
var1_rnaseq_names <- c("met20", "met21", "met22")
var2_rnaseq_names <- c("met14", "met15", "met16")

average_meth_results_directory <- "/PATH/TO/average.meth.genes.levels/mto1_vs_wt/"
DMRs_results_directory <- "/PATH/TO/Methylome.At/results/mto1_vs_wt/genome_annotation/"
RNAseq_results_directory <- "/PATH/TO/met23/mto1_vs_wt/"
main_output_directory <- "/PATH/TO/Linear_correlation/"

#-------------------------------------------------------------------#

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/multiplot_ggplot2.R")
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/yo_theme_base_ggplot2.r")
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/linear_plot_over_gene_features_DMRs_DEGs.r")

linear_plot_meth_rna(
    var1_name = var1_name,
    var2_name = var2_name,
    var1_rnaseq_names = var1_rnaseq_names,
    var2_rnaseq_names = var2_rnaseq_names,
    average_meth_results_directory = average_meth_results_directory,
    DMRs_results_directory = DMRs_results_directory,
    RNAseq_results_directory = RNAseq_results_directory,
    genes_2_keep = "filtered_by_DEGs",
    main_output_directory = main_output_directory,
    additional_plots = TRUE,
    pValues_table = TRUE
)
```
#### example output

<img
  src="https://github.com/Yo-yerush/general_scripts/blob/main/lm.stats.plot.mto1.png"
  alt="fig"
  width="100%"
/>
