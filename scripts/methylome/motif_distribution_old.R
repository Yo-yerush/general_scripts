library(org.At.tair.db)
library(rtracklayer)
library(Biostrings)

breaks = 1000
n_rdm = 1000

file_title = "test"
context = "CG"
sub_cntx =  NULL #"CCG"
ann = "Genes"
#go_type = "MF"
#go_gairORloss = "gain"
treatment = "mto1_vs_wt"
go_id_list = NULL
#go_id_list = "GO:0004123"
#go_id_list = c("GO:0090116","GO:0009664","GO:0010424")
#go_id_list = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/BSseq_results_111123/res/",treatment,"/GO_analysis/",context,"/",ann,"/",go_type,".",ann,".",context,".",go_gairORloss,".",treatment,".topGO.csv"))$GO.ID
#go_id_list = go_id_list[1:16]
DMR_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/",treatment,"/genome_annotation/",context,"/",ann,"_",context,"_genom_annotations.csv"))
#DMR_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/",treatment,"/genome_annotation/",context,"/TEG_",context,"_genom_annotations.csv"))
#tair_id_vec = NULL
tair_id_vec = unique(DMR_file$locus_tag)
#tair_id_vec = "AT4G08113"
tair_id_vec = "AT1G35480" # MRD1
#tair_id_vec = "AT3G01120" # MTO1
#tair_id_vec = "AT2G24735" # SDG21
distribution_range_values = c(NA,NA) # c(bellow,above) # from list of tairs ids, will filter all distribution less than the mention value.

#########################################################
#########################################################
# run only once
if (T) {
  #gff3 = import.gff3("/home/yoyerush/yo/TAIR10.1/GTF_file/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
  gff3 = import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
  gff3 = gff3[-which(gff3@seqnames == "NC_037304.1")]
  gff3 = gff3[-which(gff3@seqnames == "NC_000932.1")]
  
  fasta <- readDNAStringSet("C:/Users/yonye/Migal/Rachel Amir Team - General/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna")
  
  Genes <- gff3[which(gff3$gbkey == "Gene")]
  CDS <- gff3[which(gff3$gbkey == "CDS")]
  Promoters = promoters(Genes, upstream=2000, downstream=0, use.names=TRUE)
  Promoters$gbkey = "Promoter"
  
  if (ann == "Genes") {ann_type = Genes
  } else if (ann == "CDS") {
    ann_type = CDS
  } else if (ann == "Promoters") {
    ann_type = Promoters
  }
  
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
      
      if (!is.null(sub_cntx)) {
        motif_patterns = sub_cntx
        xlim = c(0,6)
        
      } else if (context == "CG") {
        motif_patterns = "CG"
        xlim = c(0,6)
        
      } else if (context == "CHG") {
        motif_patterns = c("CAG", "CTG", "CCG")
        xlim = c(0,6)
        
      } else if (context == "CHH") {
        motif_patterns = c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC")
        #"TTG", "ATG", "GTG", "TAG", "AAG", "GAG", "TGG", "AGG", "GGG")
        xlim = c(0,20)
      }
      
      #sum pattern in gene sequence
      motif_count_0 = sum(sapply(motif_patterns, function(pattern) countPattern(pattern, gene_seq[[1]])))
      # prepare distribution score
      motif_count_0 = (motif_count_0 / length(gene_seq[[1]])) * 100
      
      motif_count = append(motif_count_0,motif_count, after = F)
    })
  }
}
#########################################################
#########################################################

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
  
  # find only score above value
  if (!is.na(distribution_range_values[1]) | !is.na(distribution_range_values[2])) {
    tair2motifDis_0 = data.frame(tair = tair_id_list, motif_dis = tair_motif_count)
    tair2motifDis_0 = tair2motifDis_0[!duplicated(tair2motifDis_0$tair),]
    if (!is.na(distribution_range_values[1])) {tair2motifDis = tair2motifDis_0[tair2motifDis_0$motif_dis <= distribution_range_values,]}
    if (!is.na(distribution_range_values[2])) {tair2motifDis = tair2motifDis_0[tair2motifDis_0$motif_dis >= distribution_range_values,]}

    dis_DMR_print = DMR_file[grep(paste(unique(tair2motifDis$tair), collapse="|"), DMR_file$locus_tag),]
    tair_motif_count = unique(tair2motifDis$motif_dis)
    View(dis_DMR_print)
  }
  
  context.vis = ifelse(is.null(sub_cntx), context, sub_cntx)
#svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/",file_title,"_10kRandom",context,"_",ann,"_distribution.svg"),
#    width = 2.83, height = 3.38, family = "serif")
hist(motif_count[-1],
     main = paste0(n_rdm," random ",ann),
     xlab = paste0(context.vis," Distribution / ",ann," length (%)"),
     ylab = "Frequency",
     xlim = c(min(motif_count[-1]), max(motif_count[-1])),
     breaks = breaks)
abline(v = tair_motif_count,                      
       col = rep("red4",length(tair_motif_count)),
       lwd = 1)
#dev.off()
