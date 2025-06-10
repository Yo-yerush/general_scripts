library(dplyr)
library(rtracklayer)
library(GenomicFeatures)
#library(tidyr)
############################### Araport11
# gff3
gff3_ar11 = import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/Araport11/Jul2024/Araport11_GFF3_genes_transposons.current.gff.gz")[,1:19]
remove_seqnames = c("ChrM","ChrC")
gff3_ar11 <- gff3_ar11[!seqnames(gff3_ar11) %in% remove_seqnames]
seqlevels(gff3_ar11) = as.character(unique(seqnames(gff3_ar11)))

ar11 = gff3_ar11[,-6]

# remove '.1' and other suffix numbers from 'CDS' ids
CDS_rows = which(ar11$type == "CDS")
ar11$ID[CDS_rows] = gsub("\\.\\d","", ar11$ID[CDS_rows])

# remove duplicate 'type' and 'ID' rows
ar11$tmp = paste(as.character(seqnames(ar11)),start(ar11),end(ar11),as.character(strand(ar11)),
                 ar11$type,#ar11$ID,
                 sep = "_")
ar11 = ar11[!duplicated(ar11$tmp), -ncol(mcols(ar11))]

ar11$gene_id = ar11$ID
#ar11_genes = ar11[which(ar11$type == "gene")]

# get introns ranges
gr_introns = makeTxDbFromGRanges(ar11) %>%
  intronsByTranscript() %>% 
  unlist(.,use.names = F) %>%
  unique() %>%
  sort(., by = ~seqnames + start + end + strand)
gr_introns$type = as.factor("intron")

# Map the locus_tag based on overlap or direct ID ('locus_tag') matching
overlaps <- findOverlaps(gr_introns, ar11)
first_overlap <- tapply(seq_along(queryHits(overlaps)), queryHits(overlaps), `[`, 1)
gr_introns$gene_id <- mcols(ar11)[subjectHits(overlaps)[first_overlap], "gene_id"]
gr_introns$ID = "intron:"

# group by 'gene_id' and number the introns (considering strand direction)
gr_introns <- as.data.frame(gr_introns) %>%
  group_by(gene_id, strand) %>%
  mutate(ID = ifelse(strand == "-",
                     paste0(ID, rev(row_number())),
                     paste0(ID, row_number()))) %>%
  ungroup() %>%
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

# combine 'ar11' and 'introns'
ar11_combined = sort(c(ar11, gr_introns), by = ~seqnames + start + end + strand)
ar11_combined = ar11_combined[,c("type","gene_id","ID")]
ar11_combined$gene_id = gsub(":.*","",ar11_combined$gene_id)
ar11_combined$ID = gsub("^AT[0-5]G[0-9]{5}:","",ar11_combined$ID)


# keep only TAIR gene IDs (and specific type in 'ID' column)
ar11_combined = ar11_combined[grep("^AT[0-5]G[0-9].*[0-9]$", mcols(ar11_combined)$gene_id)]
ar11_combined$gene_id = gsub("\\.[0-9]","",ar11_combined$gene_id)



##############################################
# add gene type column
gene_type = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/Araport11/Jul2024/Araport11_gene_type",
                     sep = "\t", header = F)[-(1:4),]
names(gene_type) = c("gene_id","gene_type")


as.data.frame(ar11_combined) %>%
  merge(.,gene_type, by = "gene_id", all.x = T) %>%
  relocate(type, .before = gene_id) %>%
  rename(ID = "annotation_id") %>%
  write.csv(gzfile("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_GFF3_Araport11_Jul24.csv.gz"),
            row.names = F)
  #makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  #export.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_GFF3_Araport11_Jul24.gff3.gz")

# keep only the following annotation types
#type_vec = "^gene$|CDS|five_prime_UTR|three_prime_UTR|^pseudogene$|transposable_element_gene|lnc_RNA"
#ar11_combined = ar11_combined[grep(type_vec, mcols(ar11_combined)$type)]

export.gff3(ar11_combined,
            "C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_GFF3_Araport11_based_Jul24.gff3.gz")
import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_GFF3_Araport11_Jul24.gff3.gz")



yo = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/Methylome.At_GFF3_Araport11_Jul24.csv.gz") %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T)











# add 'intron' ID
intron_rows = which(ar11_combined$type == "intron")
for (ii in 1:length(intron_rows)) {
  i = intron_rows[ii]
  
  if (as.character(strand(ar11_combined[i])) == "+") {
    
    nearest_end = start(ar11_combined[i]) - 1
    x = ar11_combined[ar11_combined@seqnames == seqnames(ar11_combined[i]) & end(ar11_combined) == nearest_end]$ID
    ar11_combined$ID[i] = x[grep("CDS",x)][1] %>% gsub("exon:","intron:") 
    
  } else if (as.character(strand(ar11_combined[i])) == "-") {
    
    nearest_end = end(ar11_combined[i]) + 1
    x = ar11_combined[ar11_combined@seqnames == seqnames(ar11_combined[i]) & start(ar11_combined) == nearest_end]$ID
    ar11_combined$ID[i] = x[grep("CDS",x)][1] %>% gsub("exon:","intron:") 

  }
  # its take a long time.......
  cat(round((ii/length(intron_rows))*100,2),"%\n")
}






#nearest_exon_pos = ar11_combined$ID[(i-1):(i-6)] %>% grep("exon",.) %>% .[1]
#ar11_combined$ID[i] = gsub("exon:","intron:",ar11_combined$ID[i-nearest_exon_pos]) # take CDS id from row above


#nearest_exon_pos = ar11_combined$ID[(i+1):(i+6)] %>% grep("exon",.) %>% .[1]
#ar11_combined$ID[i] = gsub("exon:","intron:",ar11_combined$ID[i+nearest_exon_pos]) # take CDS id from row bellow



# if duplicate 'intron' id but different location, add suffix
ar11_combined$tmp = paste(ar11_combined$gene_id,ar11_combined$ID,sep = "_")
duplicated_introns = grep("intron",ar11_combined$tmp)


duplicated_introns = ar11_combined[intron_rows]
duplicated_introns = duplicated_introns[duplicated(duplicated_introns$tmp),]


ar11_combined[duplicated(duplicated_introns$tmp)]





# gtf
#gff3_ar11_2 = import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/Araport11/Jul2024/Araport11_GTF_genes_transposons.current.gtf.gz")

gr_introns = makeTxDbFromGRanges(gff3_ar11) %>%
  intronsByTranscript(., use.names = FALSE) %>% 
  unlist(.,use.names = F) %>%
  unique()
gr_introns$type = as.factor("intron")
# Map the locus_tag based on overlap or direct ID ('locus_tag') matching
overlaps <- findOverlaps(gr_introns, gff3_ar11)
first_overlap <- tapply(seq_along(queryHits(overlaps)), queryHits(overlaps), `[`, 1)
gr_introns$locus_tag <- mcols(gff3_ar11)[subjectHits(overlaps)[first_overlap], "locus_tag"]


mcols(gff3_ar11) <- mcols(gff3_ar11)[, c("type","ID")]
#mcols(gff3_ar11) <- mcols(gff3_ar11)[, c("type","source")]




#ar11_df = as.data.frame(gff3_ar11)
#ar11_df$tmp = paste(ar11_df$seqnames,ar11_df$start,ar11_df$end,ar11_df$strand,ar11_df$type,sep = "_")



############################### NCBI (refseq) annotation file (TAIR10.1)
gff3_refseq = import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")

gff3_refseq = renameSeqlevels(gff3_refseq, gsub("NC_003070.9","Chr1",seqlevels(gff3_refseq)))
gff3_refseq = renameSeqlevels(gff3_refseq, gsub("NC_003071.7","Chr2",seqlevels(gff3_refseq)))
gff3_refseq = renameSeqlevels(gff3_refseq, gsub("NC_003074.8","Chr3",seqlevels(gff3_refseq)))
gff3_refseq = renameSeqlevels(gff3_refseq, gsub("NC_003075.7","Chr4",seqlevels(gff3_refseq)))
gff3_refseq = renameSeqlevels(gff3_refseq, gsub("NC_003076.8","Chr5",seqlevels(gff3_refseq)))
gff3_refseq = renameSeqlevels(gff3_refseq, gsub("NC_037304.1","ChrM",seqlevels(gff3_refseq))) # MT
gff3_refseq = renameSeqlevels(gff3_refseq, gsub("NC_000932.1","ChrC",seqlevels(gff3_refseq))) # Pltd
remove_seqnames = c("ChrM","ChrC")
gff3_refseq <- gff3_refseq[!seqnames(gff3_refseq) %in% remove_seqnames]
seqlevels(gff3_refseq) = as.character(unique(seqnames(gff3_refseq)))
refseq = gff3_refseq

mcols(refseq) <- mcols(refseq)[, c("type", "locus_tag","source")]

############################### 
gr_introns = makeTxDbFromGRanges(gff3_refseq) %>%
  intronsByTranscript() %>% 
  unlist(.,use.names = F) %>%
  unique()
gr_introns$type = as.factor("intron")
# Map the locus_tag based on overlap or direct ID ('locus_tag') matching
overlaps <- findOverlaps(gr_introns, gff3_refseq)
first_overlap <- tapply(seq_along(queryHits(overlaps)), queryHits(overlaps), `[`, 1)
gr_introns$locus_tag <- mcols(gff3_refseq)[subjectHits(overlaps)[first_overlap], "locus_tag"]

############################### 
merged_gr = c(refseq,gff3_ar11,gr_introns)
merged_gr = merged_gr[order(seqnames(merged_gr), start(merged_gr), end(merged_gr),strand(merged_gr))]

merged_gr$key <- paste(as.character(seqnames(merged_gr)), 
                               start(merged_gr),
                               end(merged_gr),
                               mcols(merged_gr)$strand, 
                               mcols(merged_gr)$type, 
                               #mcols(merged_gr)$locus_tag, 
                               sep = "_")
merged_gr <- merged_gr[!duplicated(mcols(merged_gr)$key)]

type_vec = "^gene$|CDS|five_prime_UTR|three_prime_UTR|^pseudogene$|transposable_element_gene|lnc_RNA"
merged_gr = merged_gr[grep(type_vec, mcols(merged_gr)$type)]
mcols(merged_gr) = mcols(merged_gr)[,-ncol(mcols(merged_gr))]













refseq_df = as.data.frame(gff3_refseq)
refseq_df$tmp = paste(refseq_df[,1:ncol(refseq_df)],sep = "_")

refseq_df$tmp = paste(refseq_df$seqnames,refseq_df$start,refseq_df$end,refseq_df$strand,refseq_df$type,sep = "_")
refseq_df = refseq_df[!duplicated(refseq_df$tmp),]





merge.data.frame(ar11_df,refseq_df, by = "tmp", all.x = T) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T)















remove_seqnames = c("ChrM","ChrC")
gff3_refseq <- gff3_refseq[!seqnames(gff3_refseq) %in% remove_seqnames]
seqlevels(gff3_refseq) = as.character(unique(seqnames(gff3_refseq)))

gr_introns = makeTxDbFromGRanges(gff3_refseq) %>%
  intronsByTranscript() %>% 
  unlist(.,use.names = F) %>%
  unique()
gr_introns$type = as.factor("intron")

# Map the locus_tag based on overlap or direct ID ('locus_tag') matching
overlaps <- findOverlaps(gr_introns, gff3_refseq)
first_overlap <- tapply(seq_along(queryHits(overlaps)), queryHits(overlaps), `[`, 1)
gr_introns$locus_tag <- mcols(gff3_refseq)[subjectHits(overlaps)[first_overlap], "locus_tag"]

# bind and order introns into 'gff3_refseq'
gff3_refseq_final = c(gff3_refseq,gr_introns)
gff3_refseq_final = gff3_refseq_final[order(seqnames(gff3_refseq_final), start(gff3_refseq_final), end(gff3_refseq_final))]

gff3_refseq_final$key <- paste(as.character(seqnames(gff3_refseq_final)), 
                               start(gff3_refseq_final),
                               end(gff3_refseq_final),
                               mcols(gff3_refseq_final)$strand, 
                               mcols(gff3_refseq_final)$type, 
                               mcols(gff3_refseq_final)$locus_tag, 
                               sep = "_")
gff3_refseq_final <- gff3_refseq_final[!duplicated(mcols(gff3_refseq_final)$key)]
gff3_refseq_final <- gff3_refseq_final[!seqnames(gff3_refseq_final) %in% gff3_refseq_final]
mcols(gff3_refseq_final) <- mcols(gff3_refseq_final)[, c("type", "locus_tag")]




















############################### TAIR10 annotation file
gff3_tair10 = import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 gff3/TAIR10_GFF3_genes.gff")
remove_seqnames = c("ChrM","ChrC")
gff3_tair10 <- gff3_tair10[!seqnames(gff3_tair10) %in% remove_seqnames]
seqlevels(gff3_tair10) <- setdiff(seqlevels(gff3_tair10), remove_seqnames)


gr_introns = makeTxDbFromGRanges(gff3_tair10) %>%
  intronsByTranscript() %>% 
  unlist(.,use.names = F) %>%
  unique()

#gff3_tair10 = gff3_tair10[grep(type_vec, mcols(gff3_tair10)$type)]
mcols(gff3_tair10)$type = gsub("mRNA","transcript",mcols(gff3_tair10)$type)
unlist(intronsByTranscript(makeTxDbFromGRanges(gff3_tair10)))
#tair10_df = as.data.frame(gff3_tair10) %>% 
#  select(Name, Derives_from) %>%
#  filter(grepl("AT[0-5]TE.*", Derives_from)) %>%
#  rename(gene_id = Name)


overlaps <- findOverlaps(gff3_ar11, gff3_tair10)
# Create a new GRanges object based on overlaps
merged_gr <- gff3_ar11[queryHits(overlaps)]
# Add new columns from gr2 to merged_gr
mcols(merged_gr)$gene_name <- mcols(gr2)[subjectHits(overlaps), "gene_name"]











############################### merge all
merge_descriptions = merge.data.frame(ar11,uniP, by = "gene_id", all = T)
merge_descriptions = merge.data.frame(merge_descriptions,gff3_df, by = "gene_id", all = T)
merge_descriptions = merge.data.frame(merge_descriptions,tair10_df, by = "gene_id", all = T)

merge_descriptions = merge_descriptions %>%
  relocate(gene, .after = gene_id) %>%
  relocate(type, .before = Protein.families) %>% 
  relocate(Function, .before = EC.number) %>%
  relocate(db_xref, .before = note) %>%
  relocate(DOI.ID, PubMed.ID, .after = Derives_from) %>%
  select(-Organism)

write.csv(merge_descriptions,gzfile("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/AT_description_file_methylome.csv.gz"),
          row.names = F, na = "")
write.csv(merge_descriptions,"C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/AT_description_file_methylome.csv",
          row.names = F, na = "")
