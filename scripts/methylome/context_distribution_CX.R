library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(DMRcaller)
library(parallel)

for (treatment in c("mto1","mto3","35s")) {
 
  gff3_path = "/home/yoyerush/yo/TAIR10.1/GTF_file/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3"
  TE_path = "/home/yoyerush/yo/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt"
  des_path = "/home/yoyerush/yo/methylome_pipeline/Methylome.At/Methylome.At_scripts/data/AT_description_file_methylome.csv"
  import_TE_path = "/home/yoyerush/yo/methylome_pipeline/Methylome.At/Methylome.At_scripts/import_TE_file.R"
  rename_tair_path = "/home/yoyerush/yo/methylome_pipeline/Methylome.At/Methylome.At_scripts/rename_TAIR.2.TAIR.R"
  
  sample_df = read.table(paste0("/home/yoyerush/yo/methylome_pipeline/Methylome.At/samples_table/samples_table_",treatment,".txt"))
  output_df_path = "/home/yoyerush/yo/methylome_pipeline/context_distribution_CX/"
  
  
  new_path = paste0(output_df_path,treatment,"/")
  dir.create(new_path)
  
  ### create 'annotation_vec'
  {
    gff3 = import.gff3(gff3_path)
    gff3 = gff3[-which(gff3@seqnames == "NC_037304.1")]
    gff3 = gff3[-which(gff3@seqnames == "NC_000932.1")]
    
    Genes <- gff3[which(gff3$type == "gene")]
    CDS <- gff3[which(gff3$type == "CDS")]
    Promoters = promoters(Genes, upstream=2000, downstream=0, use.names=TRUE)
    Promoters$type = "Promoter"
    
    # TEs
    source(import_TE_path)
    source(rename_tair_path)
    TE_file.df = read.csv(TE_path, sep = "\t")
    TE_file.gr = import_TE(TE_file.df)
    TE_file = rename_seq_TAIR.2.TAIR(TE_file.gr, is_TAIR10.1=T)
    colnames(mcols(TE_file))[colnames(mcols(TE_file)) == "Transposon_Name"] <- "locus_tag"

    Retrotransposons <- TE_file[TE_file$Transposon_Super_Family %in% c("LTR/Copia", "LTR/Gypsy", "LINE/L1"), ]
    
    # TEG
    des_file = read.csv(des_path)
    TEG_0 = des_file[des_file$gene_model_type == "transposable_element_gene",] %>%
      select(locus_tag) %>% filter(!is.na(locus_tag))
    TEG = makeGRangesFromDataFrame(merge.data.frame(as.data.frame(Genes), TEG_0, by = "locus_tag"),
                                   keep.extra.columns = T)

    
    ######################## find over laps between gff3 and introns/UTRs
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
    
    
    annotation_vec = list(Genes = Genes,
                          Promoters = Promoters,
                          CDS = CDS,
                          Introns = Introns,
                          fiveUTRs = fiveUTRs,
                          threeUTRs = threeUTRs,
                          TEG = TEG,
                          Transposable_Elements = TE_file,
                          Retrotransposons = Retrotransposons)
  }
  
  ##################################################################
  read_function <- function(file_path) {
    readBismark(file_path)
  }
  
  # Use mclapply to read files in parallel
  read_CX <- mclapply(sample_df[,2], read_function, mc.cores = nrow(sample_df))
  gr_list <- GRangesList(read_CX)
  names(gr_list) = ave(sample_df[,1], sample_df[,1], FUN = function(x) paste0(x, ".", seq_along(x)))
  
  
  ##################################################################
  
  for (context in c("CG","CHG","CHH",
                    "CAG", "CTG", "CCG",
                    "CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC")) {
    
    # create data frame with extra 'nrow(sample_df)' columns
    ann_CX_df <- data.frame(
      type = names(annotation_vec),
      matrix(NA, nrow = length(annotation_vec), ncol = nrow(sample_df), dimnames = list(NULL, paste0("v", 1:nrow(sample_df))))
    )   
    
    names(ann_CX_df)[-1] = ave(sample_df[,1], sample_df[,1], FUN = function(x) paste0(x, ".", seq_along(x)))
    
    for (ann.l in 1:length(annotation_vec)) {
      
      ##################################################################
      ### count CX from annotations
      cx2count <- function(meth_cx) {
        
        if ((context == "CG" | context == "CHG" | context == "CHH")) {
          meth_cx = meth_cx[which(meth_cx$context == context)]
        } else {
          meth_cx = meth_cx[which(meth_cx$trinucleotide_context == context)]
        }
        meth_cx$proportion = meth_cx$readsM / meth_cx$readsN
        meth_cx = meth_cx[which(meth_cx$readsN >= 6)]
        meth_cx = meth_cx[which(meth_cx$proportion >= 0.4)]
        
        # find overlaps for all CX with annotation file (gff3)
        ## var1 (control)
        m1 <- findOverlaps(meth_cx, annotation_vec[[ann.l]])
        CX_annotation <- meth_cx[queryHits(m1)]
        return(length(CX_annotation))
      }
      
      for (var in 1:length(gr_list)) {
        ann_CX_df[ann.l,var+1] = cx2count(gr_list[[var]])
      }
    }
    
    ##################################################################
    
    v1_name = unique(sample_df[,1])[1]
    v2_name = unique(sample_df[,1])[2]
    
    v1_col = grep(v1_name, names(ann_CX_df))
    v2_col = grep(v2_name, names(ann_CX_df))
    
    for (i.row in 1:nrow(ann_CX_df)) {
      ann_CX_df$v1_m[i.row] = mean(as.numeric(ann_CX_df[i.row,v1_col]))
      ann_CX_df$v1_s[i.row] = sd(as.numeric(ann_CX_df[i.row,v1_col]))
      
      ann_CX_df$v2_m[i.row] = mean(as.numeric(ann_CX_df[i.row,v2_col]))
      ann_CX_df$v2_s[i.row] = sd(as.numeric(ann_CX_df[i.row,v2_col]))
      
    }
    ann_CX_df$ave.proportion = ann_CX_df$v2_m / ann_CX_df$v1_m
    
    names(ann_CX_df) = gsub("v1_m", paste0(v1_name,".ave"), names(ann_CX_df))
    names(ann_CX_df) = gsub("v1_s", paste0(v1_name,".sd"), names(ann_CX_df))
    names(ann_CX_df) = gsub("v2_m", paste0(v2_name,".ave"), names(ann_CX_df))
    names(ann_CX_df) = gsub("v2_s", paste0(v2_name,".sd"), names(ann_CX_df))
    
    ##################################################################
    
    write.csv(ann_CX_df, paste0(new_path,context,"_CX_count.csv"), row.names = F)
  } 
}