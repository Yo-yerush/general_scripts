library(ggplot2)
library(dplyr)
library(GenomicRanges)

source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/windowSize_for_GRanges_mcol.r")

read_mut_file <- function(mut_name) {
    x <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/mutants_figs/", mut_name, "_delta_df.csv.gz"))

    return(list(
        cg = x[x$context == "CG", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta"),
        chg = x[x$context == "CHG", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta"),
        chh = x[x$context == "CHH", ] %>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% windowSize("delta")
    ))
}

met1 <- read_mut_file("met1")
cmt2 <- read_mut_file("cmt2")
cmt3 <- read_mut_file("cmt3")
ddm1 <- read_mut_file("ddm1")

ChrPlots_CX("test_stroud_290525", list(met1$cg,cmt2$cg,cmt3$cg,ddm1$cg), c("met1","cmt2","cmt3","ddm1"), "CG", y_max = 0, y_mid = 0, y_min = -1, output_dir="C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/mutants_figs/")

ChrPlots_CX("test_stroud_290525", list(met1$chg,cmt2$chg,cmt3$chg,ddm1$chg), c("met1","cmt2","cmt3","ddm1"), "CHG", y_max = 0.2, y_mid = 0, y_min = -0.5, output_dir="C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/mutants_figs/")

ChrPlots_CX("test_stroud_290525", list(met1$chh,cmt2$chh,cmt3$chh,ddm1$chh), c("met1","cmt2","cmt3","ddm1"), "CHH", y_max = 0.05, y_mid = 0, y_min = -0.15, output_dir="C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/mutants_figs/")
