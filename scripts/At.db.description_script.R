library(dplyr)
library(org.At.tair.db)

############### 'org.At.tair.db' to data frame
group_db_fun <- function(at.df.f, delimiter=NULL, prefix=NULL, is.rev=F, new.name) {
  ## Bimap interface:
  x <- as.list(at.df.f) %>% unlist
  xx = data.frame(gene_id = names(x), db_vec = x, row.names = 1:length(x)) %>%
    filter(!is.na(db_vec)) %>%
    mutate(gene_id = ifelse(nchar(gene_id) > 9, substr(gene_id, 1, 9), gene_id))
  
  if (!is.null(prefix)) {
    xx$db_vec = paste0(prefix, xx$db_vec)
  }
  
  if (is.null(delimiter)) {
    xx = xx[!duplicated(xx$gene_id, fromLast = T), ]
    
  } else {
    xx <- xx %>%
      mutate(tmp = paste(gene_id, db_vec, sep = "_")) %>%
      distinct(tmp, .keep_all = T) %>%
      group_by(gene_id) %>%
      summarise(db_vec = if(is.rev) {
        paste(rev(db_vec), collapse = delimiter)
      } else {
        paste(db_vec, collapse = delimiter)
      }) %>%
      as.data.frame()
  }
  names(xx)[2] = new.name
  return(xx)
}
###############

symbol = group_db_fun(org.At.tairSYMBOL, NULL, NULL, is.rev=F, "Symbol")
other_symbol = group_db_fun(org.At.tairSYMBOL, " ", NULL, is.rev=T, "old_symbols")
EC = group_db_fun(org.At.tairENZYME, "; ", NULL, is.rev=F, "EC")
refseq = group_db_fun(org.At.tairREFSEQ, "; ", NULL, is.rev=F, "refseq_id")
kegg = group_db_fun(org.At.tairPATH, "; ", "ath", is.rev=F, "KEGG_pathway")
GENENAME = group_db_fun(org.At.tairGENENAME, "; ", NULL, is.rev=F, "Function")
PMID = group_db_fun(org.At.tairPMID, "; ", NULL, is.rev=F, "PMID")
AraCyc_name = group_db_fun(org.At.tairARACYCENZYME, "; ", NULL, is.rev=F, "AraCyc_name")

### 'AraCyc' is want a special treatment...
AraCyc_0 <- as.list(org.At.tairARACYC) %>% unlist
AraCyc = data.frame(gene_id = AraCyc_0, db_vec = names(AraCyc_0), row.names = 1:length(AraCyc_0)) %>%
  filter(!is.na(db_vec)) %>%
  mutate(db_vec = ifelse(grepl("\\d{1,3}$", db_vec),
                          sub("\\d{1,3}$", "", db_vec),
                         db_vec))
AraCyc = rbind(AraCyc[grep("Arabidopsis",AraCyc$db_vec),],
           AraCyc[grep("plants",AraCyc$db_vec),],
           AraCyc[-grep("Arabidopsis|plants",AraCyc$db_vec),]) %>%
  group_by(gene_id) %>%
  summarise(db_vec = paste(db_vec, collapse = "; ")) %>%
  as.data.frame()
names(AraCyc)[2] = "AraCyc"

# change the 'HTML' shit in 'AraCyc'
for (ii in c("<i>","<I>","<small>","<sub>","<sup>","<SUP>","<em>","</i>","</I>","</small>","</sub>","</sup>","</SUP>","</em>")) {
  AraCyc_name[,2] = gsub(ii, "", AraCyc_name[,2])
  AraCyc[,2] = gsub(ii, "", AraCyc[,2])
}

At.db.all = merge.data.frame(symbol, other_symbol, by = "gene_id", all = T) %>%
  merge.data.frame(., EC, by = "gene_id", all = T) %>%
  merge.data.frame(., refseq, by = "gene_id", all = T) %>%
  merge.data.frame(., kegg, by = "gene_id", all = T) %>%
  merge.data.frame(., GENENAME, by = "gene_id", all = T) %>%
  merge.data.frame(., PMID, by = "gene_id", all = T) %>%
  merge.data.frame(., AraCyc_name, by = "gene_id", all = T) %>%
  merge.data.frame(., AraCyc, by = "gene_id", all = T)



