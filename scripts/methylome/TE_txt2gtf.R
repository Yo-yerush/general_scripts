library(rtracklayer)

TE = read.csv("P:/TAIR10.1/TAIR10_Transposable_Elements.txt", sep = "\t")

TE$seqnames = TE$Transposon_Name
TE$seqnames = gsub(paste0("AT1TE.*"), "NC_003070.9",TE$seqnames)
TE$seqnames = gsub(paste0("AT2TE.*"), "NC_003071.7",TE$seqnames)
TE$seqnames = gsub(paste0("AT3TE.*"), "NC_003074.8",TE$seqnames)
TE$seqnames = gsub(paste0("AT4TE.*"), "NC_003075.7",TE$seqnames)
TE$seqnames = gsub(paste0("AT5TE.*"), "NC_003076.8",TE$seqnames)

names(TE)[2:4] = c("strand","start","end")
TE = TE[,c("seqnames","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family")]

TE$strand = gsub("true","+",TE$strand)
TE$strand = gsub("false","-",TE$strand)

TE_gr = sort(makeGRangesFromDataFrame(TE, keep.extra.columns = T))
export(TE_gr, "P:/TAIR10.1/TAIR10_TE_GTF.gtf")