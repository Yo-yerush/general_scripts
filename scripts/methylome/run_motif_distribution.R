library(org.At.tair.db)
library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(dplyr)
source("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/scripts/methylome/motif_distribution_plot.R")
source("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/scripts/methylome/motif_random_distribution.R")

#######################################
### annotated DMRs
{
  treatment = "mto1_vs_wt"
  
  n_rdm.l = 5000
  
  dir.create("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/motif_distribution/mto1/from_DMRs/",
             showWarnings = F)
  
  for (context.l in c("CG","CHG","CHH",
                      "CAG", "CTG", "CCG",
                      "CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC")) {
    
    #dir.create(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/motif_distribution/mto1/from_DMRs/",context.l),
    #           showWarnings = F)
    
    if (nchar(context.l) == 3) {
      context.ll <- context.l
      substr(context.ll, 2, 2) <- "H"
    } else {
      context.ll = context.l
    }
    
    pdf(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/motif_distribution/mto1/from_DMRs/",context.l,".pdf"),
        width = 6.75, height = 4.82, family = "serif")
    par(mfrow = c(2, 3), oma = c(0, 0, 2, 0))
    
    for (ann.l in c("Genes","Promoters","CDS","Introns","TEG","Transposable_Elements")) {
      
      DMR_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/",treatment,"/genome_annotation/",context.ll,"/",ann.l,"_",context.ll,"_genom_annotations.csv"))
      if (ann.l == "Transposable_Elements") {
        names(DMR_file) = gsub("Transposon_Name", "locus_tag", names(DMR_file))
      }
      #DMR_file = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/mto1_vs_wt/genome_annotation/CHG/TEG_CHG_genom_annotations.csv")
      try({
      random_genes = motif_random_distribution(n_rdm = n_rdm.l, ann = ann.l, context = context.l)
      
      motif_distribution(n_rdm = n_rdm.l,
                         ann = ann.l,
                         context = context.l,
                         tair_id_vec = unique(DMR_file$locus_tag),
                         main_title = ann.l,
                         motif_count = random_genes$count,
                         motif_patterns = random_genes$patterns)
      })
      cat(">")
    }
    mtext(paste0("n.random = ",n_rdm.l), outer = TRUE, cex = 1)
    dev.off()
    cat("\n")
  }
}
#######################################

#######################################
### sig. genes from 'group' excel file
{
  n_rdm.l = 5000
  
  dir.create("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/motif_distribution/mto1/from_DEG_groups/",
             showWarnings = F)
  
  for (ann.l in c("Genes","Promoters","Introns","CDS")) {
    dir.create(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/motif_distribution/mto1/from_DEG_groups/",ann.l),
               showWarnings = F)
    
    for (context.l in c("CG","CHG","CHH",
                        "CAG", "CTG", "CCG",
                        "CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC")) {
      
      try({
      random_genes = motif_random_distribution(n_rdm = n_rdm.l, ann = ann.l, context = context.l)
      })
      
      pdf(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/motif_distribution/mto1/from_DEG_groups/",ann.l,"/",context.l,"_",ann.l,".pdf"),
          width = 6.75, height = 4.82, family = "serif")
      par(mfrow = c(2, 3), oma = c(0, 0, 2, 0))
      
      for (sheet.l in c("DNA_methyltransferase","Histone_Lysine_MTs","Royal_Family_Proteins",
                        "chromatin_remodeling","RdDM_pathway","other_methylation")) {
        df_x = as.data.frame(readxl::read_xlsx(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/mto1_groups.xlsx"),
                                               sheet = sheet.l)) %>% filter(RNA_pvalue < 0.05) %>%
          select(locus_tag) %>% distinct(locus_tag, .keep_all = T)
        
        try({
        motif_distribution(n_rdm = n_rdm.l,
                           ann = ann.l,
                           context = context.l,
                           tair_id_vec = df_x$locus_tag,
                           main_title = sheet.l,
                           motif_count = random_genes$count,
                           motif_patterns = random_genes$patterns)
        })
      }
      mtext(paste0(n_rdm.l," random ",ann.l), outer = TRUE, cex = 1)
      dev.off()
      
    }
    cat("* ",ann.l," *\n")
  }
}
#######################################


#tair_id_vec = "AT4G08113"
#tair_id_vec = "AT1G35480" # MRD1
#tair_id_vec = "AT3G01120" # MTO1
#tair_id_vec = "AT2G24735" # SDG21