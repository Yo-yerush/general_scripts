
# Find all files in the specified directory
file_path <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/GSE51304_RAW"

# Get all files in the directory
files <- list.files(path = file_path, full.names = TRUE)
# Keep only files with .txt.gz suffix
files <- files[grepl("\\.txt\\.gz$", files)]

### smRNA
if (T) {
setwd("C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/Yo_files_for_IGV")
    dir.create("smRNA")
    setwd("smRNA")
    smrna_files <- files[grepl("smRNA", files)]

    for (i_file in smrna_files) {
        output_path <- gsub("txt", "wig", basename(i_file))
        txt_to_wig(i_file, output_path, smRNA = T, h3_file = NULL, bigwig = TRUE, bin_size = 50)
    }
}

### ChIP
if (T) {
    chip_files <- files[grepl("ChIP", files)]
    h3_chip_files <- chip_files[grepl("H3_ChIP", chip_files)]
    other_chip_files <- chip_files[!grepl("H3_ChIP", chip_files)]

    base_path <- "C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/GSE51304_RAW/"

    chip_files <- c(
        paste0(base_path, "GSM1242392_WT_H3_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242392_WT_H3_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242395_ddcc_H3_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242395_ddcc_H3_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242398_suvh456_H3_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242398_suvh456_H3_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242392_WT_H3_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242392_WT_H3_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1667164_WT_H3_ChIP_rep_processed.txt.gz"),
        paste0(base_path, "GSM1667164_WT_H3_ChIP_rep_processed.txt.gz"),
        paste0(base_path, "GSM1667167_ddcc_H3_ChIP_rep_processed.txt.gz"),
        paste0(base_path, "GSM1667167_ddcc_H3_ChIP_rep_processed.txt.gz"),
        paste0(base_path, "GSM1667170_suvh456_H3_ChIP_rep_processed.txt.gz"),
        paste0(base_path, "GSM1667170_suvh456_H3_ChIP_rep_processed.txt.gz")
    )
    names(chip_files) <- c(
        paste0(base_path, "GSM1242393_WT_H3K9me2_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242394_WT_H3K23ac_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242396_ddcc_H3K9me2_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242397_ddcc_H3K23ac_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242399_suvh456_H3K9me2_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242400_suvh456_H3K23ac_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242411_WT_H3K9me1_LJ_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1242412_WT_H3K9me2_LJ_ChIP_processed.txt.gz"),
        paste0(base_path, "GSM1667165_WT_H3K9me2_ChIP_rep_processed.txt.gz"),
        paste0(base_path, "GSM1667166_WT_H3K23ac_ChIP_rep_processed.txt.gz"),
        paste0(base_path, "GSM1667168_ddcc_H3K9me2_ChIP_rep_processed.txt.gz"),
        paste0(base_path, "GSM1667169_ddcc_H3K23ac_ChIP_rep_processed.txt.gz"),
        paste0(base_path, "GSM1667171_suvh456_H3K9me2_ChIP_rep_processed.txt.gz"),
        paste0(base_path, "GSM1667172_suvh456_H3K23ac_ChIP_rep_processed.txt.gz")
    )


    setwd("C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/Yo_files_for_IGV")

    dir.create("ChIP")
    setwd("ChIP")
    for (i_file in seq(length(chip_files))) {
        output_path <- gsub("txt", "wig", basename(names(chip_files)[i_file]))
        txt_to_wig(
            names(chip_files)[i_file],
            output_path,
            smRNA = F,
            h3_file = as.character(chip_files[i_file]),
            bigwig = TRUE,
            bin_size = 50
        )
    }
}



######################### methylome and RNAseq files
### mRNA seq
setwd("C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/Yo_files_for_IGV")
    dir.create("mRNA")
    setwd("mRNA")
if (T) {
    files <- list.files(path = file_path, full.names = TRUE) # mRNA
    files <- files[grepl("_mRNA", files)]
    for (i_file in files) {
        wig_data <- readr::read_lines(i_file)
        wig_data_modified <- stringr::str_replace_all(wig_data, "variableStep chrom=chr", "variableStep chrom=Chr")
        readr::write_lines(wig_data_modified, basename(i_file))
    }
}


### WGBS
setwd("C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/Yo_files_for_IGV")
    dir.create("WGBS")
    setwd("WGBS")
if (T) {
    files <- list.files(path = "C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/Jacobsen_Lab_Epigenomics_Data/GSE39901_RAW", full.names = TRUE) # 5mC
    for (i_file in files) {
        wig_data <- readr::read_lines(i_file)
        wig_data_modified <- stringr::str_replace_all(wig_data, "variableStep chrom=chr", "variableStep chrom=Chr")
        readr::write_lines(wig_data_modified, basename(i_file))
    }
}
