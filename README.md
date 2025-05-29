-----------------------------------------------------------------
# WGBS pipeline

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
### run Bismark to get just 'CX_report' file (can run without --cx option to get all ouput files)
```
.run_bismark_yo.sh -s dml3_NS/bismark_dml3_samples.txt -g /home/yoyerush/TAIR10_chr_all.fas.gz -o ./dml3_NS/bismark_results -n 32 -m 16G --cx
```

### Use Methylome.At downstream pipeline

-----------------------------------------------------------------

# Calculate 'delta' methylation levels from '.wig' files ('stroud et al. 2013' example)
#### Use [mutants compare_delta_df.r](https://github.com/Yo-yerush/general_scripts/blob/main/mutants_compare_delta_df.r) script

#### Then to create a ChrPlot use [mutants compare_delta_plot.r](https://github.com/Yo-yerush/general_scripts/blob/main/mutants_compare_delta_plot.r) script
#### ChrPlot in all context:
![fig](https://github.com/Yo-yerush/general_scripts/blob/main/ChrPlot_CG_test_stroud_290525.svg)
![fig](https://github.com/Yo-yerush/general_scripts/blob/main/ChrPlot_CHG_test_stroud_290525.svg)
![fig](https://github.com/Yo-yerush/general_scripts/blob/main/ChrPlot_CHH_test_stroud_290525.svg)

-----------------------------------------------------------------
