library(dplyr)
library(tidyr)

############################### uniProt
uniP_0 = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/uniProt_description_files/uniprotkb_arabidopsis_AND_model_organis_2024_01_25.csv")[,-1]
uniP_1 = uniP_0[grep("A",uniP_0$TAIR),]
uniP_1$TAIR = gsub(";A","<>A",uniP_1$TAIR)
uniP_1$TAIR = gsub(";","",uniP_1$TAIR)
uniP <- uniP_1 %>%
  mutate(TAIR = strsplit(as.character(TAIR), "<>")) %>%
  unnest(TAIR)
uniP = uniP[!duplicated(uniP$TAIR),] %>% 
  select(TAIR, Protein.families, everything())
names(uniP)[1] = "locus_tag"

uniP$Function..CC. = gsub("FUNCTION: ","", uniP$Function..CC.)
uniP$Tissue.specificity = gsub("TISSUE SPECIFICITY: ","", uniP$Tissue.specificity)
uniP$Induction = gsub("INDUCTION: ","", uniP$Tissue.specificity)
names(uniP) = gsub("Gene.Ontology\\.\\.","GO\\.", names(uniP))
names(uniP) = gsub("\\.\\.CC\\.","", names(uniP))

############################### Araport11
ar11 = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/Araport11_description_files/Araport11_functional_descriptions_051023_without_C_M.txt", sep = "\t")
names(ar11)[1] = "locus_tag"
ar11$locus_tag = gsub("\\.[0-9]","",ar11$locus_tag)
ar11 = ar11[!duplicated(ar11$locus_tag),]


############################### NCBI (refseq) annotation file (TAIR10.1)
gff3 = rtracklayer::import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
gff3_df = as.data.frame(gff3) %>% 
  select(locus_tag,gene,type,gene_biotype,db_xref,note) %>%
  filter(type == "gene") %>%
  filter(!grepl("Arth|DA", locus_tag)) %>%
  filter(!is.na(locus_tag))


############################### merge all
merge_descriptions = merge.data.frame(ar11,uniP, by = "locus_tag", all = T)
merge_descriptions = merge.data.frame(merge_descriptions,gff3_df, by = "locus_tag", all = T)

merge_descriptions = merge_descriptions %>%
  relocate(gene, .after = locus_tag) %>%
  relocate(type, .before = Protein.families) %>% 
  relocate(Function, .before = EC.number) %>%
  relocate(gene_model_type, db_xref, .before = note) %>%
  relocate(DOI.ID, PubMed.ID, .after = note) %>%
  select(-Organism)

write.csv(merge_descriptions,"C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/AT_description_file_methylome.csv",
          row.names = F, na = "")
