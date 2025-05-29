-----------------------------------------------------------------
-----------------------------------------------------------------

# WGBS pipeline

#### * Download fastq files from SRA and run [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/#:~:text=Bismark%20is%20a%20program%20to%20map%20bisulfite%20treated,the%20methylation%20levels%20of%20their%20samples%20straight%20away.) for *dml3* samples ([Zhejiang University](https://www.ncbi.nlm.nih.gov/sra/SRX4698864))
```
./download_fastq_from_SRA_short.sh "SRR7848067 SRR7848068 SRR7848069 SRR7848070" "/PATH/TO/dml3_NS"
```

#### * Create a sample table file (*tab* delimiter):
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
#### * Run Bismark to get only '**.CX_report.txt**' file (can run without '--cx' option to get all ouput files):
```
.run_bismark_yo.sh -s ./dml3_NS/bismark_dml3_samples.txt -g /PATH/TO/TAIR10_chr_all.fas.gz -o ./dml3_NS/bismark_results -n 32 -m 16G --cx
```

#### * Use [Methylome.At](https://github.com/Yo-yerush/Methylome.At) downstream pipeline
```

```

-----------------------------------------------------------------
-----------------------------------------------------------------

# Calculate and plot '*delta*' methylation levels from '*.wig*' files
*[Stroud et al. 2013](https://pubmed.ncbi.nlm.nih.gov/23313553/) ([SRP014726](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP014726&o=biosample_s%3Aa%3Bacc_s%3Aa))*
#### * Use [mutants compare_delta_df.r](https://github.com/Yo-yerush/general_scripts/blob/main/mutants_compare_delta_df.r) script

#### * Then to create a ChrPlot use [mutants compare_delta_plot.r](https://github.com/Yo-yerush/general_scripts/blob/main/mutants_compare_delta_plot.r) script
#### ChrPlots:
![fig](https://github.com/Yo-yerush/general_scripts/blob/main/ChrPlot_CG_test_stroud_290525.svg)
![fig](https://github.com/Yo-yerush/general_scripts/blob/main/ChrPlot_CHG_test_stroud_290525.svg)
![fig](https://github.com/Yo-yerush/general_scripts/blob/main/ChrPlot_CHH_test_stroud_290525.svg)

-----------------------------------------------------------------
-----------------------------------------------------------------

