motif_distribution <- function(n_rdm,
                               ann,
                               context,
                               #sub_cntx =  NULL, #"CCG"
                               tair_id_vec = NULL,
                               go_id_list = NULL,
                               main_title = NA,
                               distribution_range_values = "5%",#c(NA,NA) # c(bellow,above) from list of tairs ids, will filter all distribution less/above than the mention value
                               motif_count, motif_patterns
) {
  #########################################################
  breaks = 1000
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
    colnames(mcols(Transposable_Elements))[colnames(mcols(Transposable_Elements)) == "Transposon_Name"] <- "locus_tag"
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
  
  #########################################################
  # 'tair_id_list' from GO or vector of tair ids
  if (!is.null(go_id_list)) {
    # tair gene motif distribution string
    xx <- as.list(org.At.tairGO2TAIR)
    tair_id_list_b = data.frame(x=NA)
    for (i in 1:length(go_id_list)) {
      goids <- xx[go_id_list[i]]
      tair_id_list_0 = data.frame(x = goids[[1]])
      tair_id_list_b = rbind(tair_id_list_b,tair_id_list_0)
    }
    tair_id_list_v = data.frame(x = unique(tair_id_list_b$x[-1]))
    
    tair_id = data.frame(x = unique(DMR_file$locus_tag))
    
    tair_id_list = unique(merge.data.frame(tair_id_list_v,tair_id, by = "x")[,1])
  }
  
  if (!is.null(tair_id_vec)) {tair_id_list = tair_id_vec}
  
  #########################################################
  # motif count in selected genes
  #message("motif count in ",length(tair_id_list)," selected genes...")
  tair_motif_count = NA
  for (tair_no in 1:length(tair_id_list)) {
    
    ann_tair_gr = ann_type[which(ann_type$locus_tag == tair_id_list[tair_no]),]
    chr_name = unique(as.character(ann_tair_gr@seqnames))
    chr_pos = grep(chr_name, names(fasta))
    
    start_pos = min(ann_tair_gr@ranges@start)
    end_pos = max(start_pos + ann_tair_gr@ranges@width - 1)
    gene_seq = subseq(fasta[chr_pos,], start = start_pos, end = end_pos)
    
    if (ann_tair_gr@strand@values == "-") {
      gene_seq =  reverseComplement(gene_seq)
    }
    
    #sum pattern in gene sequence
    motif_count_0 = sum(sapply(motif_patterns, function(pattern) countPattern(pattern, gene_seq[[1]])))
    # prepare distribution score
    motif_count_0 = (motif_count_0 / length(gene_seq[[1]])) * 100
    
    tair_motif_count = append(motif_count_0,tair_motif_count, after = F)
  }
  tair_motif_count = tair_motif_count[-1]
  
  #########################################################
  # find only score above value
  if (distribution_range_values == "5%") {
    lower_cutoff <- quantile(motif_count, 0.05)
    upper_cutoff <- quantile(motif_count, 0.95)
    tair_5perx_list = (tair_motif_count <= lower_cutoff | tair_motif_count >= upper_cutoff)
    tair_motif_count_5perc = tair_motif_count[tair_5perx_list]
    #print(tair_id_list[tair_5perx_list])
    
  } else {
  if (!is.na(distribution_range_values[1]) | !is.na(distribution_range_values[2])) {
    tair2motifDis_0 = data.frame(tair = tair_id_list, motif_dis = tair_motif_count)
    tair2motifDis_0 = tair2motifDis_0[!duplicated(tair2motifDis_0$tair),]
    if (!is.na(distribution_range_values[1])) {tair2motifDis = tair2motifDis_0[tair2motifDis_0$motif_dis <= distribution_range_values,]}
    if (!is.na(distribution_range_values[2])) {tair2motifDis = tair2motifDis_0[tair2motifDis_0$motif_dis >= distribution_range_values,]}
}
    dis_DMR_print = DMR_file[grep(paste(unique(tair2motifDis$tair), collapse="|"), DMR_file$locus_tag),]
    tair_motif_count = unique(tair2motifDis$motif_dis)
    View(dis_DMR_print)
  }
  
  #########################################################
  new_main_title = ifelse(is.na(main_title), paste0(n_rdm," random ",ann),
                          gsub("_"," ",main_title))
  # context.vis = ifelse(is.null(sub_cntx), context, sub_cntx)
  #svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/",file_title,"_10kRandom",context,"_",ann,"_distribution.svg"),
  #    width = 2.83, height = 3.38, family = "serif")
  hist_data = hist(motif_count,
                   main = new_main_title,
                   xlab = paste0(context," Distribution / ",ann," length (%)"),
                   ylab = "Frequency",
                   xlim = c(min(motif_count), max(motif_count)),
                   breaks = breaks)
  if (distribution_range_values == "5%") {
  abline(v = tair_motif_count,                      
         col = rep("#7d9ec950",length(tair_motif_count)),
         #col = rainbow(length(tair_motif_count)),
         lwd = 1)
  abline(v = tair_motif_count_5perc,                      
         col = rep("#94101070",length(tair_motif_count_5perc)),
         lwd = 1)
  lines(hist_data$mids, hist_data$counts, type = "h", lwd = 1, col = "gray10")
  } else {
    abline(v = tair_motif_count,                      
           col = rep("#94101080",length(tair_motif_count)),
           lwd = 1)
  }
  #dev.off()
}