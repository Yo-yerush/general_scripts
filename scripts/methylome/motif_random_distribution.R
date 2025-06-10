motif_random_distribution <- function(n_rdm,
                               ann,
                               context) {
  # gff3 annotations
  gff3 = import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
  gff3 = gff3[-which(gff3@seqnames == "NC_037304.1")]
  gff3 = gff3[-which(gff3@seqnames == "NC_000932.1")]
  
  Genes <- gff3[which(gff3$type == "gene")]
  CDS <- gff3[which(gff3$type == "CDS")]
  Promoters = promoters(Genes, upstream=2000, downstream=0, use.names=TRUE)
  Promoters$type = "Promoter"
  
  #if (ann == "Transposable_Elements") {
    source("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/scripts/methylome/import_TE_file.R")
    source("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/scripts/methylome/rename_TAIR.2.TAIR.R")
    TE_file.df = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt", sep = "\t")
    TE_file.gr = import_TE(TE_file.df)
    Transposable_Elements = rename_seq_TAIR.2.TAIR(TE_file.gr, is_TAIR10.1=T)
  #}
  
  #if (ann == "TEG") {
    des_file = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/AT_description_file_methylome.csv")
    TEG_0 = des_file[des_file$gene_model_type == "transposable_element_gene",] %>%
      select(locus_tag) %>% filter(!is.na(locus_tag))
    TEG = makeGRangesFromDataFrame(merge.data.frame(as.data.frame(Genes), TEG_0, by = "locus_tag"),
                                   keep.extra.columns = T)
  #}
  
  ######################## find over laps between gff3 and introns/UTRs
  if (ann == "Introns" | ann == "fiveUTRs" | ann == "threeUTRs") {
    gff3_TxDb = makeTxDbFromGRanges(gff3)
    
    Introns_0 = unlist(intronsByTranscript(gff3_TxDb))[,!1:3]
    # find overlaps for Introns with annotation file (gff3)
    m <- findOverlaps(Introns_0, gff3)
    Introns <- Introns_0[queryHits(m)]
    mcols(Introns) <- cbind.data.frame(
      mcols(Introns),
      mcols(gff3[subjectHits(m)]))
    Introns = unique(Introns)
    Introns$type = "intron"
    
    fiveUTRs_0 = unlist(fiveUTRsByTranscript(gff3_TxDb))[,!1:3]
    # find overlaps for fiveUTRs with annotation file (gff3)
    m <- findOverlaps(fiveUTRs_0, gff3)
    fiveUTRs <- fiveUTRs_0[queryHits(m)]
    mcols(fiveUTRs) <- cbind.data.frame(
      mcols(fiveUTRs),
      mcols(gff3[subjectHits(m)]))
    fiveUTRs = unique(fiveUTRs)
    fiveUTRs$type = "five_prime_UTR"
    
    threeUTRs_0 = unlist(threeUTRsByTranscript(gff3_TxDb))[,!1:3]
    # find overlaps for threeUTRs with annotation file (gff3)
    m <- findOverlaps(threeUTRs_0, gff3)
    threeUTRs <- threeUTRs_0[queryHits(m)]
    mcols(threeUTRs) <- cbind.data.frame(
      mcols(threeUTRs),
      mcols(gff3[subjectHits(m)]))
    threeUTRs = unique(threeUTRs)
    threeUTRs$type = "three_prime_UTR"
    
    
    ann_type = switch(ann,
                      "Introns" = Introns,
                      "fiveUTRs" = fiveUTRs,
                      "threeUTRs" = threeUTRs
    )
  } else {
    ann_type = switch(ann,
                      "Genes" = Genes,
                      "CDS" = CDS,
                      "Promoters" = Promoters,
                      "TEG" = TEG,
                      "Transposable_Elements" = Transposable_Elements
    )
  }
  
  
  fasta <- readDNAStringSet("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna")
  
  ######################################################### 
  # motif count in random genes
  #message("motif count in ",n_rdm," random genes...")
  motif_count = NA
  for (gene_no in sample(1:length(ann_type), n_rdm)) {
    try({
      chr_name = as.character(ann_type[gene_no,]@seqnames)
      chr_pos = grep(chr_name, names(fasta))
      
      start_pos = ann_type[gene_no,]@ranges@start
      end_pos = start_pos + ann_type[gene_no,]@ranges@width - 1
      gene_seq = subseq(fasta[chr_pos,], start = start_pos, end = end_pos)
      #matchPattern("AAA",gene_seq[[1]])
      if (ann_type[gene_no,]@strand@values == "-") {
        gene_seq =  reverseComplement(gene_seq)
      }
      
      if (context == "CHG") {
        motif_patterns = c("CAG", "CTG", "CCG")
        
      } else if (context == "CHH") {
        motif_patterns = c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC")
        #"TTG", "ATG", "GTG", "TAG", "AAG", "GAG", "TGG", "AGG", "GGG")
      } else {
        motif_patterns = context
      }
      
      #sum pattern in gene sequence
      motif_count_0 = sum(sapply(motif_patterns, function(pattern) countPattern(pattern, gene_seq[[1]])))
      # prepare distribution score
      motif_count_0 = (motif_count_0 / length(gene_seq[[1]])) * 100
      
      motif_count = append(motif_count_0,motif_count, after = F)
    })
  }
  motif_count = motif_count[-1]
  
  return(list(count = motif_count,
              patterns = motif_patterns))
}