library(org.At.tair.db)
library(rtracklayer)
library(Biostrings)
if (F) {
  dir.create("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/motif_distribution/random_genes")
  dir.create("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/motif_distribution/DMRs_list")
  dir.create("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/motif_distribution/GO_genes")
}

#gff3 = import.gff3("/home/yoyerush/yo/TAIR10.1/GTF_file/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
gff3 = import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
gff3 = gff3[-which(gff3@seqnames == "NC_037304.1")]
gff3 = gff3[-which(gff3@seqnames == "NC_000932.1")]

fasta <- readDNAStringSet("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna")

Genes <- gff3[which(gff3$gbkey == "Gene")]
CDS <- gff3[which(gff3$gbkey == "CDS")]
Promoters = promoters(Genes, upstream=2000, downstream=0, use.names=TRUE)
Promoters$gbkey = "Promoter"

######################################################################
######################################################################
###### for random genes
random_genes = function(context,n_rdm,ann,breaks) {
  if (ann == "Genes") {
    ann_type = Genes
  } else if (ann == "CDS") {
    ann_type = CDS
  } else if (ann == "Prmoters") {
    ann_type = Promoters
  }
  try({
    motif_count = NA
    for (gene_no in sample(1:length(ann_type), n_rdm)) {
      try({
        chr_name = as.character(ann_type[gene_no,]@seqnames)
        chr_pos = grep(chr_name, names(fasta))
        
        start_pos = ann_type[gene_no,]@ranges@start
        end_pos = start_pos + ann_type[gene_no,]@ranges@width - 1
        gene_seq = subseq(fasta[chr_pos,], start = start_pos, end = end_pos)
        #matchPattern("AAA",gene_seq[[1]])
        
        if (context == "CG") {
          motif_patterns = "CG"
          xlim = c(0,8)
          
        } else if (context == "CHG") {
          motif_patterns = c("CAG", "CTG", "CCG", "CGG") # CGG because its CCG in reverse string
          xlim = c(0,8)
          
        } else if (context == "CHH") {
          motif_patterns = c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC",
                             "TTG", "ATG", "GTG", "TAG", "AAG", "GAG", "TGG", "AGG", "GGG")
          xlim = c(15,40)
        }
        
        #sum pattern in gene sequence
        motif_count_0 = sum(sapply(motif_patterns, function(pattern) countPattern(pattern, gene_seq[[1]])))
        # prepare distribution score
        motif_count_0 = (motif_count_0 / length(gene_seq[[1]])) * 100
        
        motif_count = append(motif_count_0,motif_count, after = F)
      })
    }
    svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/BSseq_results_100823/motif_distribution/random_genes/","random_",context,"_",ann,"_distribution.svg"),
        width = 2.83, height = 3.38, family = "serif")
    hist(motif_count[-1],
         main = paste0(n_rdm," random ",ann),
         xlab = paste0(context," Distribution / ",ann," length"),
         ylab = "Frequency",
         xlim = xlim,
         breaks = breaks)
    abline(v = mean(motif_count[-1]),                      
           col = "red4",
           lwd = 3)
    dev.off()
    #motif_mean = mean(motif_count[-1])
  })
  try({dev.off()})
}

for (context_loop in c("CG","CHG","CHH")) {
  random_genes(context_loop,10000,"Genes", breaks = 1000)
  random_genes(context_loop,10000,"CDS", breaks = 50)
  random_genes(context_loop,10000,"Prmoters", breaks = 1000)
  
  message(paste0("*\n*\n*\ndone *",context_loop,"* context in *random_genes* function\n*\n*\n*"))
}

######################################################################
######################################################################
###### for DMRs
DMRs_genes = function(treatment,context,ann,breaks) {
  if (ann == "Genes") {
    ann_type = Genes
  } else if (ann == "CDS") {
    ann_type = CDS
  } else if (ann == "Prmoters") {
    ann_type = Promoters
  }
  try({
    DMR_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/BSseq_results_100823/",treatment,"/genome_annotation/",context,"/",ann,"_",context,"_genom_annotations.csv"))
    tair_id = unique(DMR_file$locus_tag)
    
    tair2score = data.frame(tair = tair_id, score = NA) # to get score of genes
    motif_count = NA
    for (tair_no in 1:length(tair_id)) {
      try({
        chr_name = unique(as.character(ann_type[which(ann_type$locus_tag == tair_id[tair_no])]@seqnames))
        chr_pos = grep(chr_name, names(fasta))
        
        start_pos = min(ann_type[which(ann_type$locus_tag == tair_id[tair_no]),]@ranges@start)
        end_pos = max(start_pos + ann_type[which(ann_type$locus_tag == tair_id[tair_no]),]@ranges@width - 1)
        gene_seq = subseq(fasta[chr_pos,], start = start_pos, end = end_pos)
        #matchPattern("AAA",gene_seq[[1]])
        
        if (context == "CG") {
          motif_patterns = "CG"
          xlim = c(0,8)
          
        } else if (context == "CHG") {
          motif_patterns = c("CAG", "CTG", "CCG", "CGG") # CGG because its CCG in reverse string
          xlim = c(0,8)
          
        } else if (context == "CHH") {
          motif_patterns = c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC",
                             "TTG", "ATG", "GTG", "TAG", "AAG", "GAG", "TGG", "AGG", "GGG")
          xlim = c(15,40)
        }
        
        #sum pattern in gene sequence
        motif_count_0 = sum(sapply(motif_patterns, function(pattern) countPattern(pattern, gene_seq[[1]])))
        # prepare distribution score
        motif_count_0 = (motif_count_0 / length(gene_seq[[1]])) * 100
        tair2score$score[tair_no] = motif_count_0 # to get score of genes
        
        motif_count = append(motif_count_0,motif_count, after = F)
      })
    }
    svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/BSseq_results_100823/motif_distribution/DMRs_list/",treatment,"_",context,"_",ann,"_distribution.svg"),
        width = 2.83, height = 3.38, family = "serif")
    hist(motif_count,
         main = paste0(gsub("_"," ",treatment)," - ",ann),
         xlab = paste0(context," Distribution / ",ann," length"),
         ylab = "Frequency",
         xlim = xlim,
         breaks = breaks)
    abline(v = mean(motif_count[-1]),                      
           col = "red4",
           lwd = 3)
    dev.off()
    #motif_mean = mean(motif_count[-1])
  })
  try({dev.off()})
}

for (treatment_loop in c("mto1_vs_wt_2","mto3_vs_wt_2")) {
  for (context_loop in c("CG","CHG","CHH")) {
    DMRs_genes(treatment_loop,context_loop,"Genes", breaks = 150)
    DMRs_genes(treatment_loop,context_loop,"CDS", breaks = 100)
    DMRs_genes(treatment_loop,context_loop,"Prmoters", breaks = 150)
  }
  message(paste0("*\n*\n*\ndone *",context_loop,"* context in *DMRs_genes* function\n*\n*\n*"))
}

# to test specific genes
if (F) {
  treatment = "mto1_vs_wt_2"
  context = "CHG"
  ann = "Genes"
  breaks = 150
  
  # run to <tair2score>
  genes_of_interest = data.frame(a = tair2score[tair2score$score >= 6.25,1])
  gene_des = read.table("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/gene_description_20140101_TAIR10.txt", sep = "\t", fill = T)
  names(gene_des)[1] = "a"
  gene_des$a = gsub("\\.[0-9]","",gene_des$a)
  gene_des_df = merge.data.frame(genes_of_interest, gene_des, by = "a")
  
}
######################################################################
######################################################################
###### for specific gene list (TAIR ID)
list_genes = function(treatment,context,ann,go_id_list,expr) {
  
  if (ann == "Genes") {
    ann_type = Genes
  } else if (ann == "CDS") {
    ann_type = CDS
  } else if (ann == "Prmoters") {
    ann_type = Promoters
  }
  
  try({
    xx <- as.list(org.At.tairGO2TAIR)
    tair_id_list_b = data.frame(x=NA)
    for (i in 1:length(go_id_list)) {
      goids <- xx[go_id_list[i]]
      tair_id_list_0 = data.frame(x = goids[[1]])
      tair_id_list_b = rbind(tair_id_list_b,tair_id_list_0)
    }
    tair_id_list_v = data.frame(x = unique(tair_id_list_b$x[-1]))
    
    DMR_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/BSseq_results_100823/",treatment,"/genome_annotation/",context,"/",ann,"_",context,"_genom_annotations.csv"))
    tair_id = data.frame(x = unique(DMR_file$locus_tag))
    
    tair_id_list = unique(merge.data.frame(tair_id_list_v,tair_id, by = "x")[,1])
    
    #tair_id_list = tair_id[grep(paste(tair_id_list_v, collapse = "|"), tair_id)]
  })
  
  try({
    motif_count = NA
    for (tair_no in 1:length(tair_id_list)) {
      try({
        chr_name = unique(as.character(ann_type[which(ann_type$locus_tag == tair_id_list[tair_no])]@seqnames))
        chr_pos = grep(chr_name, names(fasta))
        
        start_pos = min(ann_type[which(ann_type$locus_tag == tair_id_list[tair_no]),]@ranges@start)
        end_pos = max(start_pos + ann_type[which(ann_type$locus_tag == tair_id_list[tair_no]),]@ranges@width - 1)
        gene_seq = subseq(fasta[chr_pos,], start = start_pos, end = end_pos)
        #matchPattern("AAA",gene_seq[[1]])
        
        if (context == "CG") {
          motif_patterns = "CG"
          xlim = c(0,8)
          
        } else if (context == "CHG") {
          motif_patterns = c("CAG", "CTG", "CCG", "CGG") # CGG because its CCG in reverse string
          xlim = c(0,8)
          
        } else if (context == "CHH") {
          motif_patterns = c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC",
                             "TTG", "ATG", "GTG", "TAG", "AAG", "GAG", "TGG", "AGG", "GGG")
          xlim = c(15,40)
        }
        
        #sum pattern in gene sequence
        motif_count_0 = sum(sapply(motif_patterns, function(pattern) countPattern(pattern, gene_seq[[1]])))
        # prepare distribution score
        motif_count_0 = (motif_count_0 / length(gene_seq[[1]])) * 100
        
        motif_count = append(motif_count_0,motif_count, after = F)
      })
    }
    svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/BSseq_results_100823/motif_distribution/GO_genes/","GO_",expr,"_",treatment,"_",context,"_",ann,"_distribution.svg"),
        width = 2.83, height = 3.38, family = "serif")
    hist(motif_count,
         main = paste0(gsub("_"," ",treatment)," - ",ann),
         xlab = paste0(context," Distribution / ",ann," length"),
         ylab = "Frequency",
         xlim = xlim,
         breaks = length(motif_count[-1]))
    abline(v = mean(motif_count[-1]),                      
           col = "red4",
           lwd = 3)
    dev.off()
    #motif_mean = mean(motif_count[-1])
  })
  try({dev.off()})
}

for (treatment_loop in c("mto1_vs_wt_2","mto3_vs_wt_2")) {
  for (context_loop in c("CG","CHG","CHH")) {
    for (ann_loop in c("Genes","CDS","Prmoters")) {
      for (go_procces_loop in c("BP","CC","MF")) {
        for (gainORloss_lopp in c("gain","loss")) {
          
          list_genes(treatment_loop, context_loop, ann_loop,
                     read.csv(paste0(
                       "C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/BSseq_results_100823/GO_analysis/",treatment_loop,"/",context_loop,"/",ann_loop,"/",treatment_loop,".",go_procces_loop,".",ann_loop,".",context_loop,".",gainORloss_lopp,".topGO.csv")
                     )$GO.ID,
                     paste(go_procces_loop,gainORloss_lopp, sep = "_"))
          
          message(paste0("*\n*DMRs_genes* function done: ",
                         paste(treatment_loop,context_loop,ann_loop,go_procces_loop,gainORloss_lopp, sep = "-")
                         ))
          
        }
      }
    }
  }
}

if (F) {
  treatment = "mto1_vs_wt_2"
  context = "CHG"
  ann = "Genes"
  go_procces = "BP"
  gainORloss = "gain"
  go_id_list = read.csv(paste0(
    "C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/BSseq_results_100823/GO_analysis/",treatment,"/",context,"/",ann,"/",treatment,".",go_procces,".",ann,".",context,".",gainORloss,".topGO.csv")
  )$GO.ID
  expr = paste("score_6",go_procces,gainORloss, sep = "_")
}
######################################################################
######################################################################
###### for specific gene (TAIR ID)
if (F) {
  
  tair_id = "AT1G71920"
  
  dmr_gene = ann_type[which(ann_type$locus_tag == tair_id)]
  chr_name = as.character(dmr_gene[,]@seqnames)
  chr_pos = grep(chr_name, names(fasta))
  
  start_pos = dmr_gene[,]@ranges@start
  end_pos = start_pos + dmr_gene[,]@ranges@width - 1
  gene_seq = subseq(fasta[chr_pos,], start = start_pos, end = end_pos)
  
  dmr_motif_count_0 = (
    # + string
    countPattern("CAG",gene_seq[[1]]) +
      countPattern("CTG",gene_seq[[1]]) +
      countPattern("CCG",gene_seq[[1]]) +
      # - string
      countPattern("GTC",gene_seq[[1]]) +
      countPattern("GAC",gene_seq[[1]]) +
      countPattern("GGC",gene_seq[[1]])
  ) / length(gene_seq[[1]]) * 100
  dmr_motif_count_0
}

