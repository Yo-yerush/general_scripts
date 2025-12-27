Workflows:

1. Download SRA files (as '*.fastq*')
2. WGBS pipeline - using [Bismark software](https://www.bioinformatics.babraham.ac.uk/projects/bismark/#:~:text=Bismark%20is%20a%20program%20to%20map%20bisulfite%20treated,the%20methylation%20levels%20of%20their%20samples%20straight%20away.)
3. [Methylome.At](https://github.com/Yo-yerush/Methylome.At) downstream pipeline

4. Plot **delta** methylation levels from '*.wig*' or '*.CX_report.txt*' files
5. genePlots
6. long/short TEs
7. sub-context analysis - add it

8. RNAseq pipeline - using [RSEM software](https://github.com/deweylab/RSEM)
9. RNAseq downstream pipeline - fix it
10. DEGs - GO term summary (Term and its offspring)
11. DMRs-DEGs correlations
12. Create 'expression + DMRs' integrated tables from WGBS and RNAseq data
13. Classify the integrated table into gene groups (epigenetic, metabolism, stress-related, etc.)
14. Methylated TEs near DEGs - add it

15. KEGG pathways - enrichment and viewPlots - add it
16. GO analysis - add it
17. REVIGO scripts - add it
    
18. Metabolome Box-Plots for grouped or replicates data (include normalization)

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
``` bash
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
<img
  src="https://github.com/Yo-yerush/general_scripts/blob/main/ChrPlot_test_stroud_all.svg"
  alt="fig"
  width="100%"
/>

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

## TE size and distance from centromer

```r
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/TEs_size_n_distance.R")

samples_table <- read.table("/home/yoyerush/yo/methylome_pipeline/Methylome.At_180825/Methylome.At/samples_table/samples_table_mto1.txt")
TE_context_list <- TE_delta_meth(samples_table, n.cores = 6)

#############################################
## TE methylation levels (delta) and size
TE_delta_list <- list(
    te_size_plot(TE_context_list, "CG"),
    te_size_plot(TE_context_list, "CHG"),
    te_size_plot(TE_context_list, "CHH")
)

ggsave(
    filename = paste0(output_path, "TE_size_delta_scatter.png"),
    plot = gridExtra::grid.arrange(grobs = TE_delta_list, nrow = 1, ncol = 3),
    width = 10500,
    height = 2500,
    units = "px",
    dpi = 1200
)

#############################################
## TE methylation levels (delta) and distance from centromer
TE_distance <- distance_from_centromer(TE_context_list, TE_gr, window_size = 1e6)

ggsave(
    filename = paste0(output_path, "centromere_distance_methylations.png"),
    plot = TE_distance$plot,
    width = 3500,
    height = 2500,
    units = "px",
    dpi = 1200
)

```

distance from centromer data frame:

> head(TE_distance$df)
```
     te_id distance      avg_meth context width superfamily
AT1TE00010 14833064  0.0000000000     CHG    79   LTR/Copia
AT1TE00010 14833064 -0.0005698006     CHH    79   LTR/Copia
AT1TE00010 14833064  0.0000000000      CG    79   LTR/Copia
AT1TE00020 14828054  0.0000000000     CHG   126 RC/Helitron
AT1TE00020 14828054  0.0317174054     CHH   126 RC/Helitron
AT1TE00020 14828054  0.0000000000      CG   126 RC/Helitron
```

#### examples output

<div style="display: flex; gap: 10px;">
  <img
    src="https://github.com/Yo-yerush/general_scripts/blob/main/TE_size_delta_scatter.png"
    alt="fig"
    width="74%"
  />
  <img
    src="https://github.com/Yo-yerush/general_scripts/blob/main/centromere_distance_methylations.png"
    alt="fig"
    width="24%"
  />
</div>

-----------------------------------------------------------------
-----------------------------------------------------------------

## long/short TEs

```r
#### use long/short TEs script from mto1 paper
#### create new delta_metaplots script
delta_metaplots(list_type = "all", superfamily = "", max_value = 0.05, context_legend = T)
delta_metaplots("long", context_legend = T)
delta_metaplots("short")
delta_metaplots("superfamilies", "Gypsy")
```

-----------------------------------------------------------------
-----------------------------------------------------------------

## RNAseq pipeline 
#### Run [RSEM](https://github.com/deweylab/RSEM) for *dml3* samples ([Zhejiang University](https://www.ncbi.nlm.nih.gov/sra/SRX4698864))

 Download *Arabidopsis* reference genome ([TAIR10](https://www.arabidopsis.org/))
 ```bash
 cd /PATH/TO
 wget -O TAIR10_chr_all.fas.gz https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
 wget -O Araport11_GTF_genes_transposons.current.gtf.gz https://www.arabidopsis.org/api/download-files/download?filePath=Genes/Araport11_genome_release/Araport11_GTF_genes_transposons.20241001.gtf.gz

 gunzip TAIR10_chr_all.fas.gz
 gunzip Araport11_GTF_genes_transposons.current.gtf.gz
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

Run RSEM to get only '**.genes.results**' file
* *run without '--genes_results' option to get all output files*
```bash
wget https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/run_rsem_yo.sh
chmod +x run_rsem_yo.sh
./run_rsem_yo.sh -s samples_table.txt -g TAIR10_chr_all.fa.gz -a Araport11_GTF_genes_transposons.current.gtf.gz -n 32 --genes_results
```
#### Use downstream pipeline script (utilyzing [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) package)
```r
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/DEseq_fun_yo.R")
path <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23"
description_file <- "Methylome.At_description_file.csv.gz"

deseq_fc("dml3_vs_wt", "At", c(1:3,7:9), "dml3", "wt", path = path, description_file = description_file)

#### fix path and add colData argument
```
-----------------------------------------------------------------
-----------------------------------------------------------------

## DEGs - GO term summary (Term and its offspring)
After runing the RNAseq downstream pipeline, load the 'all genes' results file and choose GO ID(s) to summary
* *will summary the GO term and its offstpring terms*
* *each GO ID output one row, which include the Term, count of total unique genes related to this term (and its offspring); significants, and up-/down- regulated; and the precentage of significats compare to total genes*
```r
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/methylome/GO_DEGs_offspring_summary.R")

rnaseq_res <- read.csv("PATH/TO/all_genes_results_mto1_vs_wt.csv")
GO_offspring_summary(rnaseq_res, c("GO:0006952","GO:0006950","GO:0009607","GO:0009628"))
```
This will output:
```
Parent_GO_ID                   Category Total Upregulated Downregulated Significant Percentage
GO:0006952             defense response  1268         229           294         523    41.25 %
GO:0006950           response to stress  3155         593           738        1331    42.19 %
GO:0009607  response to biotic stimulus  1146         231           258         489    42.67 %
GO:0009628 response to abiotic stimulus  1984         385           539         924    46.57 %
```
-----------------------------------------------------------------
-----------------------------------------------------------------

## DEGs - GO term abiotic stress groups -  enrichment test
After runing the RNAseq downstream pipeline, load the 'all genes' results file
* *will enrich (fisher) the abiotic GO terms and its offstpring terms*
* use *dataset = "up"* or *dataset = "down"* argument for *up-/down- regulated DEGs*
```r
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/methylome/GO_abiotic_stress_enrichment_test.R")

rnaseq_res <- read.csv("PATH/TO/all_genes_results_mto1_vs_wt.csv")
abiotic_stress_enrichment_test(rnaseq_res, sub_title_prefix = "mto1 vs. wt")
```
This will plot:

<img
  src="https://github.com/Yo-yerush/general_scripts/blob/main/GO_stress_response_enrichmrnt_results_plot_mto1_all_DEGs.svg"
  alt="fig"
  width="65%"
/>

```r
abiotic_stress_enrichment_test(rnaseq_res, print_as_table = TRUE)
```
This will output:
```
         Test_name Sig_in_term Total_in_term Fold_enrichment  P_value Odds_ratio Significant
      Cold stress          42           118            1.40 8.47e-03       1.63        TRUE
   Osmotic stress         163           440            1.46  2.9e-08       1.75        TRUE
      Salt stress           1             2            1.97    0.443       2.95       FALSE
Water deprivation          40           102            1.55  1.4e-03       1.91        TRUE
       DNA damage          46           303            0.60    1.000       0.52       FALSE
 Oxidative stress          70           169            1.63 3.33e-06       2.10        TRUE
      UV-B stress           1             6            0.66    0.827       0.59       FALSE
         Wounding          12            36            1.32    0.180       1.47       FALSE
      Heat stress          24            88            1.08    0.378       1.10       FALSE
```
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


-----------------------------------------------------------------
-----------------------------------------------------------------

## Create 'expression + DMRs' integrated tables from WGBS and RNAseq results
*by default, return **non-unique** IDs data frame if there is **>1 DMR** within gene-body/promoter*
*if 'make_it_unique = TRUE', so the DMR's log2FC column will returned with non-numeric values*
* *DMRs files outputed from 'Methylome.At' pipeline with 'genome_annotations' directory*
* *the RNAseq result file must have gene_id+log2FC+padj+pvalue column (1-4 positions; will change to this names)*

```r
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/expression_n_DMRs_merge_results.R")

RNAseq_results <- read.csv("/PATH/TO/all_genes_results_mut_vs_wt.csv") # can filter data frame

merged_df <- expression_n_DMRs_merge_results(
    DMRs_results_directory = "/PATH/TO/Methylome.At/results/mut_vs_wt/genome_annotation/",
    RNAseq_results = RNAseq_results,
    make_it_unique = FALSE,
    designed_excel_output = NULL, # write output prefix to produce excel file
    DMR_sep = ", ", # if unique
    empty_DMR = "" # if unique
)

write.csv(merged_df, "/PATH/TO/merge_rnaseq_methylome_tables_mto1_vs_wt.csv", row.names = F)
```

-----------------------------------------------------------------
-----------------------------------------------------------------

## Classify the integrated table into gene groups (epigenetic, metabolism, stress-related, etc.)
*return **unique** IDs data frame if there is **>1 DMR** within annotation*

```r
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/expression_n_DMRs_merged_into_groups.R")
expression_n_DMRs_merged_into_groups(
  treatment = "mto1_vs_wt",
  merged_results_df = merged_df, # the merged data framte from the section above
  save_excel = TRUE,
  save_volcano_plots = TRUE,
  save_pie_plots = FALSE,
  combine_volcano_plots = TRUE, # if false, will output seperate plots for each group
  make_it_unique = FALSE, # if TRUE, each gene can obtain in one group only (the first by order)
  datasets_dir = "https://raw.githubusercontent.com/Yo-yerush/RA_lab_db/refs/heads/main/Arabidopsis/Groups_genes_list",
  output_dir = "./"
)
```
-----------------------------------------------------------------
-----------------------------------------------------------------

## GCMS plots
Will take tables of values and output an boxplot, to replicates or grouped samples (using '*group_lines = TRUE*' argument)
```r
source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/metabolome_boxplot_yo.R")

GCMS_results_df <- read.csv("PATH/TO/mto1_GCMS.csv")

GCMS_results_df[, 1] <- tools::toTitleCase(GCMS_results_df[, 1])
uni_compounds <- unique(GCMS_results_df[, 1])

all_plots <- list()
for (i in 1:length(uni_compounds)) {
    all_plots[[i]] <- metabolome_boxplot(
        comp.name = uni_compounds[i],
        df = GCMS_results_df,
        exp = "mto1",
        ctrl = "wt",
        log_norm = TRUE,
        group_lines = FALSE,
        p_as_star = TRUE,
        x_cex = 12,
        box_col = ifelse(
            (uni_compounds[i] == "Methionine" | uni_compounds[i] == "TOTAL AA"),
            "#629ecf", "gray60"
        ),
        path_for_save_file = NULL
    )
}

ggsave(
    paste0("PATH/TO/plots/GCMS_mto1_vs_wt_replicates_", Sys.Date(), ".svg"),
    plot = gridExtra::grid.arrange(grobs = all_plots, nrow = 3, ncol = 6),
    width = 11,
    height = 7
)
```
#### example output
<img
  src="https://github.com/Yo-yerush/general_scripts/blob/main/GCMS_mto1_vs_wt_rep_2025-12-24.svg"
  alt="fig"
  width="100%"
/>