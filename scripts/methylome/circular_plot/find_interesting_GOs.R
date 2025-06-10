library(dplyr)
library(writexl)
library(stringr)
library(RColorBrewer)

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

all_interesting_tairs = do.call(rbind, all_DFs)[,-1]
row.names(all_interesting_tairs) = 1:nrow(all_interesting_tairs)

############ legend ############
svg(paste0("P:/yonatan/methionine/circular_plot_res/corr_TE_exp/terms_legend.svg"), width = 4, height = 6, family = "serif")
plot.new()
legend("top",
       legend=names(all_DFs),
       col=brewer.pal(n=length(all_DFs)+1, "Set1")[-6],
       pch=26, lty = 1, lwd = ,
       bty="n",
       title="Terms")

dev.off()