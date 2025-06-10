library(dplyr)
library(tidyr)


# edit uniProt description file
uniP_0 = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/uniProt_description_files/uniprotkb_arabidopsis_AND_model_organis_2024_01_25.csv")[,-1]
uniP_1 = uniP_0[grep("A",uniP_0$TAIR),]
uniP_1$TAIR = gsub(";A","<>A",uniP_1$TAIR)
uniP_1$TAIR = gsub(";","",uniP_1$TAIR)
uniP <- uniP_1 %>%
  mutate(TAIR = strsplit(as.character(TAIR), "<>")) %>%
  unnest(TAIR)
uniP = uniP[!duplicated(uniP$TAIR),] %>% select(TAIR, Protein.families, everything())
names(uniP)[1] = "locus_tag"

# Araport11
ar11 = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/Araport11_description_files/Araport11_functional_descriptions_051023_without_C_M.txt", sep = "\t")
names(ar11)[1] = "locus_tag"
ar11$locus_tag = gsub("\\.[0-9]","",ar11$locus_tag)
ar11 = ar11[!duplicated(ar11$locus_tag),]


#tair2refseq = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10/TAIR10 NCBI mapping files/TAIR10_NCBI_REFSEQ_mapping_RNA", sep = "\t", header = F)[,-1]
#names(tair2refseq) = c("transcript_id", "locus_tag")
#tair2refseq$locus_tag = gsub("\\.[0-9]","",tair2refseq$locus_tag)

gff3 = rtracklayer::import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
gff3_filtered = data.frame(gff3$type,
                           #transcript_id = gsub("transcript:","",gff3$ID),
                           locus_tag = gff3$locus_tag,
                           gene = gff3$gene)
gff3_filtered = gff3_filtered[gff3_filtered$gff3.type == "gene", -1] %>% distinct(locus_tag, .keep_all = T)

#merge_descriptions = merge.data.frame(tair2refseq,ar11, by = "locus_tag", all.x = T)
merge_descriptions = merge.data.frame(gff3_filtered,ar11, by = "locus_tag", all.x = T)
merge_descriptions = merge.data.frame(merge_descriptions,uniP, by = "locus_tag", all = T)

write.csv(merge_descriptions,"C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/AT_description_file_230624.csv", row.names = F)








if (F) {
  uniP = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/uniProt_description_files/uniprotkb_arabidopsis_AND_model_organis_2024_01_25.csv")[,-1]
  names(uniP)[1] = "locus_tag"
  uniP = uniP[grep("AT", uniP$locus_tag),]
  uniP$locus_tag = gsub(";$", "",uniP$locus_tag)
  uniP_split <- uniP %>%
    separate_rows(locus_tag, sep = ";")
  uniP_split = uniP_split[!duplicated(uniP_split$locus_tag),]
}