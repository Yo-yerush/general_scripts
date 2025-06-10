library(dplyr)
library(org.At.tair.db)
library(topGO)
rnaseq_res = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv")
rnaseq_res_sig = rnaseq_res %>% filter(padj < 0.05)
#############################################################
#############################################################
#############################################################
# op 1
x <- org.At.tairGO
# Get the TAIR gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

tair2GO = data.frame(locus_tag = unique(names(xx)),
                     go_term = NA)

for (i in 1:nrow(tair2GO)) {
  tair2GO$go_term[i] = paste(unique(names(xx[[i]])), collapse = ", ")
}
write.table(tair2GO,"remove.tair2go.txt", sep = "\t", col.names = F, row.names = F, quote = F)

geneID2GO <- readMappings(file = "remove.tair2go.txt") 
file.remove("remove.tair2go.txt")

# Assuming you have a vector of interesting gene IDs
interestingGenes <- rnaseq_res_sig$locus_tag
# Universe of genes: all genes that could have been identified in your experiment
allGenes <- rnaseq_res$locus_tag # Replace with all gene IDs from your experiment

# Create a named vector indicating if a gene is interesting (1) or not (0)
geneList <- factor(as.integer(allGenes %in% interestingGenes))
names(geneList) <- allGenes

# Define the function for gene selection
geneSelectionFun <- function(allGenes) {
  return(allGenes)
}

myGOdata <- new("topGOdata", ontology="BP", allGenes=geneList, geneSel = function(x) x == 1, 
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher")

scores <- score(resultFisher)
specificGO <- "GO:0009086"  # Replace with your specific GO term
scores[specificGO]


allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher")

myterms = c("GO:0032259", "GO:0009086", "GO:0071555")
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms)) {
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  mygenesforterm <- paste(mygenesforterm, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm))
}


#############################################################
#############################################################
#############################################################
# op 2
fisher_results <- function(term_name, RNA = RNA_df) {
  
  term_name_0 = term_name
  term_name = gsub("\\[","\\\\\\[", term_name)
  term_name = gsub("\\]","\\\\\\]", term_name)
  term_df = rbind(RNA[grep(term_name, RNA$Gene.Ontology..biological.process.),],
                  RNA[grep(term_name, RNA$Gene.Ontology..molecular.function.),],
                  RNA[grep(term_name, RNA$Gene.Ontology..cellular.component.),],
                  RNA[grep(term_name, RNA$Protein.families),])
  
  n.upregulated = term_df %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0) %>% nrow()
  n.downregulated = term_df %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0) %>% nrow()
    
  a = term_df %>% filter(padj < 0.05) %>% nrow() # significant in group
  #b = term_df %>% filter(padj >= 0.05) %>% nrow() # *not* significant in group
  b = nrow(term_df) - a # *not* significant in group
  c = (RNA %>% filter(padj < 0.05) %>% nrow()) - a # total significant (without group)
  d = nrow(RNA) - c - a - b # total *not* significant (without group)
  
  
  contingency_table = matrix(c(a,c,b,d),
                             nrow = 2,
                             dimnames = list(c("term_df", "Not term_df"),
                                             c("significant", "Not significant")))
  
  fisher = fisher.test(contingency_table, alternative = "greater")
  enrichment_pValue = fisher$p.value
  enrichment_score = as.numeric(fisher$estimate)
  #enrichment_score = (a/b)/(c/d)
  
  return(data.frame(term = term_name_0, score = enrichment_score, pValue = enrichment_pValue,
                    upregulated = n.upregulated, downregulated = n.downregulated))
}

View(
  rbind(
    fisher_results("DNA methylation [GO:0006306]"),
    fisher_results("histone H3-K9 methylation [GO:0051567]"),
    fisher_results("DNA methylation-dependent heterochromatin formation [GO:0006346]"),
    fisher_results("methionine biosynthetic process [GO:0009086]"),
    fisher_results("S-adenosylmethionine biosynthetic process [GO:0006556]"),
    fisher_results("S-adenosylmethionine-dependent methyltransferase activity [GO:0008757]"),
    fisher_results("pyridoxal phosphate binding [GO:0030170]"), # mto
    fisher_results("methylation [GO:0032259]"),
    fisher_results("chromatin [GO:0000785]"),
    fisher_results("chromatin organization [GO:0006325]"),
    fisher_results("chromatin remodeling [GO:0006338]"),
    fisher_results("SAM-binding methyltransferase superfamily"),
    fisher_results("Class I-like SAM-binding methyltransferase superfamily"),
    fisher_results("Class V-like SAM-binding methyltransferase superfamily"),
    fisher_results("stress"),
    fisher_results("cell wall organization [GO:0071555]"),
    fisher_results("epigenetic regulation of gene expression [GO:0040029]"),
    fisher_results("Histone-lysine methyltransferase family"),
    fisher_results("RNA methylation [GO:0001510]"),
    fisher_results("secondary metabolic process [GO:0019748]"),
    fisher_results("flavonoid biosynthetic process [GO:0009813]"),
    fisher_results("methylglyoxal catabolic process to D-lactate via S-lactoyl-glutathione [GO:0019243]"),
    fisher_results("pectin biosynthetic process [GO:0045489]"),
    fisher_results("pectin"),
    fisher_results("DNA methylation on cytosine within a CG sequence [GO:0010424]"),
    fisher_results("DNA methylation")
  ) %>% arrange(pValue)
)
