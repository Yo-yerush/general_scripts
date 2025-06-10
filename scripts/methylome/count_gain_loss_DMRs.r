
res = list()
for (i.feature in c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs")) {

    CG = read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/", i.feature, "_CG_genom_annotations.csv"))

    CHG = read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/", i.feature, "_CHG_genom_annotations.csv"))

    CHH = read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/", i.feature, "_CHH_genom_annotations.csv"))

    total <- rbind(
        read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/", i.feature, "_CG_genom_annotations.csv")),
        read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/", i.feature, "_CHG_genom_annotations.csv")),
        read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/", i.feature, "_CHH_genom_annotations.csv"))
    )

    cg_gain = nrow(CG[CG$regionType == "gain", ])
    cg_total = nrow(CG[CG$regionType, ])

    chg_gain = nrow(CHG[CHG$regionType == "gain", ])
    chg_total = nrow(CHG[CHG$regionType, ])

    chh_gain = nrow(CHH[CHH$regionType == "gain", ])
    chh_total = nrow(CHH[CHH$regionType, ])

    total_gain = nrow(total[total$regionType == "gain", ])
    total_total = nrow(total[total$regionType, ])

    res[[i.feature]] = data.frame(
        x = c("CG", "CHG", "CHH", "total"),
        gain = round(c(
            (cg_gain / cg_total) * 100,
            (chg_gain / chg_total) * 100,
            (chh_gain / chh_total) * 100,
            (total_gain / total_total) * 100
        ), 1)
    )
}

res
