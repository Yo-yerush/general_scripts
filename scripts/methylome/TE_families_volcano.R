library(GenomicRanges)

#####################
### volcano its not good!!!!
### dont do it


TE = read.csv("P:/TAIR10.1/TAIR10_Transposable_Elements.txt", sep = "\t")


TE_DMRs_overlap <- function(TE_super_family, TE_df=TE, is.tairs=F) {

  TE_df$Transposon_ID = TE_df$Transposon_Name
  TE_df$Transposon_Name = gsub(paste0("AT1TE.*"), "NC_003070.9",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT2TE.*"), "NC_003071.7",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT3TE.*"), "NC_003074.8",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT4TE.*"), "NC_003075.7",TE_df$Transposon_Name)
  TE_df$Transposon_Name = gsub(paste0("AT5TE.*"), "NC_003076.8",TE_df$Transposon_Name)
  
  names(TE_df)[1:4] = c("seqnames","strand","start","end")
  TE_df$strand = ifelse(TE_df$strand == "true","+","-")
  TE_df$strand = "*"
  TE_gr = makeGRangesFromDataFrame(TE_df, keep.extra.columns = T)

  TE_gr_super_family = TE_gr[TE_gr$Transposon_Super_Family == TE_super_family]
  
  
  CG = read.csv("P:/yonatan/methionine/methylome_23/BSseq_results_191223/mto1_vs_wt/DMRs_CG_mto1_vs_wt.csv")[,c("seqnames","start","end","width","strand","direction")]
  CHG = read.csv("P:/yonatan/methionine/methylome_23/BSseq_results_191223/mto1_vs_wt/DMRs_CHG_mto1_vs_wt.csv")[,c("seqnames","start","end","width","strand","direction")]
  CHH = read.csv("P:/yonatan/methionine/methylome_23/BSseq_results_191223/mto1_vs_wt/DMRs_CHH_mto1_vs_wt.csv")[,c("seqnames","start","end","width","strand","direction")]
  
  DMRs = rbind(CG,CHG,CHH)
  DMRs$tmp_col = paste(DMRs$seqnames,DMRs$start,DMRs$end, sep = "XXX")
  DMRs = DMRs[!duplicated(DMRs$tmp_col),]
  DMRs = DMRs[,-grep("tmp_col", names(DMRs))]
  DMRs_gr = makeGRangesFromDataFrame(DMRs, keep.extra.columns = T)
  
  # DMRs overlap with super family
  super_family_TE_overlaps <- findOverlaps(TE_gr_super_family, DMRs_gr)
  super_family_TE_overlap_ranges <- DMRs_gr[subjectHits(super_family_TE_overlaps)]
  mcols(super_family_TE_overlap_ranges) <- cbind(mcols(super_family_TE_overlap_ranges),mcols(TE_gr_super_family)[queryHits(super_family_TE_overlaps),])
  super_family_TE_overlap_ranges = super_family_TE_overlap_ranges[!duplicated(super_family_TE_overlap_ranges$Transposon_ID)]
  
  volcano_df = data.frame(Transposon_ID = super_family_TE_overlap_ranges$Transposon_ID,
                          Transposon_ID = super_family_TE_overlap_ranges$Transposon_Family,
                          Transposon_ID = super_family_TE_overlap_ranges$direction)
  return(list(fisher = fisher.test(contingency_table)$p.value,
              overlaped = num_overlap_with_DMRs_in_family,
              annotated = num_family,
              contingency_table = contingency_table,
              direction_family = direction_family_df,
              direction_super_family = direction_super_family_df))
}
