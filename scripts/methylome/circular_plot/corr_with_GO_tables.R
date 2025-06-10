library(dplyr)
library(writexl)
library(stringr)
library(RColorBrewer)
library(rtracklayer)
library(circlize)

# part1 = find interesting GOs
{
  yo=read.csv("P:/yonatan/methionine/rnaseq_23/description_file_161123.csv")
  yoGO = yo[,grep("locus_tag|Gene.Ontology",names(yo))]
  
  if (F) {
    GO2Tair = rbind(data.frame(locus_tag = yoGO$locus_tag ,term = yoGO$Gene.Ontology..biological.process.),
                    data.frame(locus_tag = yoGO$locus_tag ,term = yoGO$Gene.Ontology..molecular.function.),
                    data.frame(locus_tag = yoGO$locus_tag ,term = yoGO$Gene.Ontology..cellular.component.))
    GO2Tair = GO2Tair[!GO2Tair$term == "",] %>% na.omit()
    GO2Tair = as.data.frame(do.call(rbind, apply(GO2Tair, 1, function(x) {
      do.call(expand.grid, strsplit(x, "; "))
    })))
    GO2Tair$locus_tag = as.character(GO2Tair$locus_tag)
    GO2Tair$term = as.character(GO2Tair$term)
    GO2Tair$tmp = paste(GO2Tair$locus_tag, GO2Tair$term, sep = " XXX ")
    GO2Tair = GO2Tair[!duplicated(GO2Tair$tmp),-3]
    GO2Tair$term = str_extract(GO2Tair$term, "GO:\\d+")
    write.csv(GO2Tair, "P:/TAIR10.1/GO_2_Tair.csv", row.names = F)
  }
  GO2Tair = read.csv("P:/TAIR10.1/GO_2_Tair.csv")
  
  
  bp_df = data.frame(x=yoGO$Gene.Ontology..biological.process.) %>% na.omit()
  bp_df = as.data.frame(do.call(rbind, apply(bp_df, 1, function(x) {
    do.call(expand.grid, strsplit(x, "; "))
  })))
  bp = as.character(unique(bp_df$x))
  
  mf_df = data.frame(x=yoGO$Gene.Ontology..molecular.function.) %>% na.omit()
  mf_df = as.data.frame(do.call(rbind, apply(mf_df, 1, function(x) {
    do.call(expand.grid, strsplit(x, "; "))
  })))
  mf = as.character(unique(mf_df$x))
  
  cc_df = data.frame(x=yoGO$Gene.Ontology..cellular.component.) %>% na.omit()
  cc_df = as.data.frame(do.call(rbind, apply(cc_df, 1, function(x) {
    do.call(expand.grid, strsplit(x, "; "))
  })))
  cc = as.character(unique(cc_df$x))
  
  GO_df = rbind(data.frame(type = rep("biological_process",length(bp)),term = bp),
                data.frame(type = rep("molecular_function",length(mf)),term = mf),
                data.frame(type = rep("cellular_component",length(cc)),term = cc))
  
  all_methylation = GO_df[grep("methylation",GO_df$term, ignore.case = TRUE),]
  all_methyltransferase = GO_df[grep("methyltransferase",GO_df$term, ignore.case = TRUE),]
  all_histone = GO_df[grep("histone",GO_df$term, ignore.case = TRUE),]
  all_DNAc5 = GO_df[grep("DNA \\(cytosine-5", GO_df$term, ignore.case = TRUE),]
  all_chromatin = GO_df[grep("chromatin",GO_df$term, ignore.case = TRUE),]
  all_CellWall = GO_df[grep("cell wall",GO_df$term, ignore.case = TRUE),]
  all_methionine = GO_df[grep("methionine",GO_df$term, ignore.case = TRUE),]
  all_others = GO_df[grep("epigenetic|methylated|methyl-Cp|COMPASS",GO_df$term, ignore.case = TRUE),]
  
  all_DFs = list(methyltransferase = all_methyltransferase,
                 DNA_5mC = all_DNAc5,
                 histone = all_histone,
                 chromatin = all_chromatin,
                 methylation = all_methylation,
                 cell_wall = all_CellWall,
                 methionine = all_methionine,
                 others = all_others)
  #write_xlsx(all_DFs, "P:/yonatan/methionine/interesting_GO_terms/interesting_GO_terms.xlsx")
  
  # extract the go ids for the list
  for (i in 1:length(all_DFs)) {
    all_DFs[[i]]$term = str_extract(all_DFs[[i]]$term, "GO:\\d+")
    for (ii in 1:nrow(all_DFs[[i]])) {
      all_DFs[[i]]$locus_tag[ii] = paste(GO2Tair[grep(all_DFs[[i]]$term[ii], GO2Tair$term), "locus_tag"],
                                         collapse = "; ")
    }
    all_DFs[[i]] = all_DFs[[i]][,-1]
    all_DFs[[i]] = as.data.frame(do.call(rbind, apply(all_DFs[[i]], 1, function(x) {
      do.call(expand.grid, strsplit(x, "; "))
    })))
    all_DFs[[i]]$term = as.character(all_DFs[[i]]$term)
    all_DFs[[i]]$locus_tag = as.character(all_DFs[[i]]$locus_tag)
    all_DFs[[i]]$col = brewer.pal(n=length(all_DFs)+1, "Set1")[-6][i]
  }
  
  
  # to save all_df data frame with description
  if (T) {
    all_DFs_des_0 = list()
    for (i in 1:length(all_DFs)) {
      all_DFs_des_0[[i]] = data.frame(locus_tag = unique(all_DFs[[i]][, "locus_tag"]))
      all_DFs_des_0[[i]] = merge.data.frame(all_DFs_des_0[[i]], yo, by = "locus_tag")
      all_DFs_des_0[[i]] = all_DFs_des_0[[i]][, -grep("locus_tag", names(all_DFs_des_0[[i]]))]
    }
    names(all_DFs_des_0) = names(all_DFs)
  }
  
  all_interesting_tairs = do.call(rbind, all_DFs)[,-1]
  row.names(all_interesting_tairs) = 1:nrow(all_interesting_tairs)
}


# part2 - circular corr plot
{
  trimm_Chr <- function(gr_obj) {
    remove_seqnames = c("NC_000932.1","NC_037304.1")
    gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
    seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
    return(sort(gr_obj))
  }
  
  ############# gff3 fNULL############# gff3 file #############
  gff3 = import.gff3("P:/TAIR10.1/GCF_000001735.4_TAIR10.1_genomic_trimmed.gff3")
  gff3.trimmed = trimm_Chr(gff3)
  #genes = gff3.trimmed[which(gff3.trimmed$type == "gene")]
  transcripts = gff3.trimmed[which(gff3.trimmed$type == "transcript")]
  
  ############# TE file #############
  TE = read.csv("P:/TAIR10.1/TAIR10_Transposable_Elements.txt", sep = "\t")
  TE_4_dens = TE[,c(1,3,4)]
  TE_4_dens$Transposon_Name = gsub(paste0("AT1TE.*"), "NC_003070.9",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT2TE.*"), "NC_003071.7",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT3TE.*"), "NC_003074.8",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT4TE.*"), "NC_003075.7",TE_4_dens$Transposon_Name)
  TE_4_dens$Transposon_Name = gsub(paste0("AT5TE.*"), "NC_003076.8",TE_4_dens$Transposon_Name)
  
  
  for (treatment in c("mto1","mto3")) {
    for (ann in c("genes","promoters")) {
      
      all_DFs_des_list = list(CG=NA, CHG=NA, CHH=NA)
      for (context in c("CG","CHG","CHH")) {
        
        ann.m = ifelse(ann == "genes", "Genes","Promoters")
        ############# methylome and RNAseq file filtered by corr files #############
        corr_file = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/corr_with_methylations/",treatment,"/",context,"/",ann,".corr.",context,".",treatment,".csv"))
        corr_file = corr_file[corr_file$pval < 0.05, c("transcript_id","locus_tag","cor")]
        #corr_file = corr_file[(corr_file$cor>0.8 | corr_file$cor<(-0.8)), c("transcript_id","locus_tag")]
        
        RNA_file = read.csv(paste0("P:/yonatan/methionine/rnaseq_23/met23/",treatment,"_vs_wt/all.transcripts.",treatment,"_vs_wt.DE.csv"))
        #RNA_file = RNA_file[RNA_file$padj < 0.05, c("transcript_id","log2FoldChange")]# %>% na.omit()
        RNA_file = RNA_file[, c("transcript_id","log2FoldChange","padj")]
        names(RNA_file)[2:3] = c("RNA_log2FC","RNA_padj")
        
        transcript2merge = as.data.frame(transcripts)[,c("seqnames","start","end","ID")]
        names(transcript2merge)[4] = "transcript_id"
        transcript2merge$transcript_id = gsub("transcript:","",transcript2merge$transcript_id)
        RNA_transcriptsLoc = merge.data.frame(transcript2merge,RNA_file, by = "transcript_id")
        
        
        RNA_filtered = merge.data.frame(RNA_file, corr_file, by = "transcript_id")
        #RNA_filtered$RNA_padj[RNA_filtered$RNA_padj > 0.05 | is.na(RNA_filtered$RNA_padj)] = "nf"
        RNA_filtered = RNA_filtered[RNA_filtered$RNA_padj < 0.05,]
        
        meth_file = read.csv(paste0("P:/yonatan/methionine/methylome_23/BSseq_results_191223/",treatment,"_vs_wt/genome_annotation/",context,"/",ann.m,"_",context,"_genom_annotations.csv"))
        meth_file$log2FoldChange = log2(meth_file$proportion2 / meth_file$proportion1)
        meth_file = meth_file[,c("locus_tag","seqnames","start","end","log2FoldChange","pValue")]
        names(meth_file)[5:6] = c("meth_log2FC","meth_padj")
        
        Chr_df = merge.data.frame(meth_file, RNA_filtered, by = "locus_tag")#[,c("seqnames","start","end","meth_log2FC","RNA_log2FC","locus_tag")]
        
        all_DFs_des = all_DFs_des_0
        for (i in 1:length(all_DFs_des)) {
          all_DFs_des[[i]] = merge.data.frame(Chr_df, all_DFs_des[[i]], by = "transcript_id")
          all_DFs_des[[i]]$tmp = paste(all_DFs_des[[i]]$locus_tag,all_DFs_des[[i]]$seqnames,all_DFs_des[[i]]$start,all_DFs_des[[i]]$end, sep = "XX")
          all_DFs_des[[i]] = all_DFs_des[[i]][!duplicated(all_DFs_des[[i]]$tmp), -grep("tmp", names(all_DFs_des[[i]]))]
          if (nrow(all_DFs_des[[i]]) != 0) {
            all_DFs_des[[i]]$Term = names(all_DFs_des)[i]
            all_DFs_des[[i]] = all_DFs_des[[i]] %>% select(Term, everything())
            all_DFs_des[[i]] = all_DFs_des[[i]] %>% relocate(cor, .before = grep("meth_log2FC",names(all_DFs_des[[i]])))
            #all_DFs_des[[i]] = all_DFs_des[[i]][,c(ncol(all_DFs_des[[i]]), 1:ncol(all_DFs_des[[i]])-1)]
          } else {
            all_DFs_des[[i]] = NA 
          }
        }
        all_DFs_des = all_DFs_des[!is.na(all_DFs_des)]
        all_DFs_des = do.call(rbind, all_DFs_des)
        
        all_DFs_des_list[[context]] = all_DFs_des
      }
      write_xlsx(all_DFs_des_list, paste0("P:/yonatan/methionine/circular_plot_res/corr_TE_exp/",treatment,"_corr_",ann,"_circular_plots.xlsx"))
    }
  }
}