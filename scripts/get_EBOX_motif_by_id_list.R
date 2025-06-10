#source("P:/yonatan/scripts/genom_fasta_promoter_2k_EBOX.R")
#source("P:/yonatan/scripts/genom_fasta_promoter_3k_EBOX.R")
peel_file = read.csv("P:/yonatan/RNAseq_yonatan_2021/DESeq2/peel/statistics/all_high_VS_low/peel.all.all_high_VS_low.DE.csv")
id_list = sample(peel_file$transcript_id,1000)
#id_list = new_table_p.val$transcript_id
#log2FC_list = new_table_p.val$log2FoldChange


i=1
while (i <= length(id_list)) {
  try(
    genom_fasta_by_refseq(id_list[i])
  )
  i=i+1
}






## upload E-BOX motif
dir = "random_genes_large_281122"
path = paste0("P:/yonatan/saved.fasta.R/EBOX_motif_peel_rnaseq_sig_genes/",dir,"/")
files = list.files(path)

if (file.exists(paste0(path,"EBOX_motifs.txt")) == T) {files = files[-grep("EBOX_motifs.txt",files)]}
id = gsub(" genom_seq_promoter_EBOX.*","", files)
df = data.frame(id = gsub(" genom_seq_promoter_EBOX.*","", files),
                bhlh_motif = "")#,
#                up_or_down = ""
                #log2FC = new_table_p.val$log2FoldChange,
                #padj = new_table_p.val$padj
#                )
i=1
while (i <= length(files)) {
  new_path = paste0(path,files[i],"/")
  
  if (file.exists(paste0(new_path,list.files(new_path, pattern = "EBOX_motif_class"))) == T) {
    
    motif_file = read.table(paste0(new_path,list.files(new_path, pattern = "EBOX_motif_class")),
                            header = T)[,1]
  } else {motif_file = "no_motif"}

  df[i,2] = paste0(motif_file, collapse=";")

  if (exists("log2FC_list") == T) {
  if (log2FC_list[i] > 0) {
    df$up_or_down[i] = "+"
  } else {df$up_or_down[i] = "-"}
  }
  
  i=i+1
}

#df[order(df$up_or_down,decreasing = T),]
#write.table(df[order(df$up_or_down,decreasing = T),], paste0(path,"EBOX_motifs.txt"),
            #row.names = F, quote = F, sep = "\t")
write.table(df, paste0(path,"EBOX_motifs.txt"),
            row.names = F, quote = F, sep = "\t")



#############################

