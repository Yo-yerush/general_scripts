sdh_GO = function(feature.name,
                  up_down, # "up", "down",  or "all"
                  n.nodes = 6) {
  
  ori_path = getwd()
  
  path = "P:/yonatan/RNAseq_yonatan_2021/SDH_variants_RNAseq_yonatan/SDH_gene_puttern_hunter/" # path for output files
  setwd(path)
  
  dir.create(paste0(path,"/GO_analysis/",feature.name))
  path.to.save.files = paste0(path,"/GO_analysis/",feature.name)
  
#  dir.create(paste0(path,"/GO_analysis/",feature.name,"/up"))
#  dir.create(paste0(path,"/GO_analysis/",feature.name,"/down"))
  
library(topGO)
library(Rgraphviz)
library(dplyr)

#treatment = "bHLH94_Vs_EV"

uniprot_pome = read.csv("P:/yonatan/RNAseq_yonatan_2021/DESeq2/uniprot/uniprot-pomegranate.csv")[,c(1,3)]
uniprot_2_pro = read.csv("P:/yonatan/RNAseq_yonatan_2021/DESeq2/uniprot/make.description.and.GO/uni.by.XP.csv")[,-c(2:4)]
pro_2_GO = merge.data.frame(uniprot_2_pro, uniprot_pome, all.x = T)[,-1]
pro_2_GO = pro_2_GO[!pro_2_GO$Gene.ontology.IDs == "",]
pro_2_GO$Gene.ontology.IDs = gsub("; ",",",pro_2_GO$Gene.ontology.IDs)
pro_2_GO$protein_id = noquote(pro_2_GO$protein_id)

trnscript_GO_id = read.csv(paste0("P:/yonatan/RNAseq_yonatan_2021/SDH_variants_RNAseq_yonatan/SDH_gene_puttern_hunter/",feature.name,"/",feature.name,"_correlation_description_table.csv"))

if (up_down == "up") {
  trnscript_GO_id = trnscript_GO_id[trnscript_GO_id$correlation >= 0.7,]
} else if (up_down == "down") {
  trnscript_GO_id = trnscript_GO_id[trnscript_GO_id$correlation <= -0.7,]
} else if (up_down == "all") {
trnscript_GO_id = rbind(trnscript_GO_id[trnscript_GO_id$correlation >= 0.7,],
                        trnscript_GO_id[trnscript_GO_id$correlation <= -0.7,])
}

#trnscript_GO_id = na.omit(trnscript_GO_id)

diff_exp_pro = trnscript_GO_id[,7]


write.table(pro_2_GO, paste0(path.to.save.files,"/remove_1.txt"), 
            row.names = F, col.names = F,quote = 2, sep = "\t")
write.table(diff_exp_pro, paste0(path.to.save.files,"/remove_2.txt"),
            row.names = F, col.names = F, quote = F)
##### μϊχο!!!!!!!!
if (file.info(paste0(path.to.save.files,"/remove_1.txt"))$size == 0) {
  file.remove(c(paste0(path.to.save.files,"/remove_1.txt"),
                paste0(path.to.save.files,"/remove_2.txt")))
  stop()
}
#####
geneID2GO <- readMappings(file = paste0(path.to.save.files,"/remove_1.txt"))   
geneUniverse <- names(geneID2GO)   
genesOfInteres <- read.table(file = paste0(path.to.save.files,"/remove_2.txt"), header = F)
genesOfInteres <- as.character(genesOfInteres$V1)

geneList<- factor(as.integer(geneUniverse %in% genesOfInteres))
#geneList<- as.integer(factor(geneUniverse %in% genesOfInteres))
names(geneList) <- geneUniverse


file.remove(c(paste0(path.to.save.files,"/remove_1.txt"),
              paste0(path.to.save.files,"/remove_2.txt")))


#### BP Ontology  ####
myGOdataBP <- new("topGOdata", ontology= "BP", allGenes=geneList, 
                  annot= annFUN.gene2GO, gene2GO= geneID2GO, nodeSize = 5)

sgBP <- sigGenes(myGOdataBP)
str(sgBP)
numSigGenes(myGOdataBP)

resultFisherBP <- runTest(myGOdataBP, algorithm = "weight01", statistic = "fisher")
allResBP <- GenTable(myGOdataBP, weightFisher= resultFisherBP, 
                     orderBy= "weightFisher", ranksOf= "weightFisher", topNodes=length(resultFisherBP@score))
allResBP$weightFisher = as.numeric(allResBP$weightFisher)
allResBP$Term = gsub(",", ";", allResBP$Term)

allResBP_view = allResBP[allResBP$weightFisher <= 0.05,]

write.csv(allResBP_view, paste0(path.to.save.files,"/",feature.name,"_corr_BP_",up_down,"_topGO.csv"), quote = F, row.names = F)

#showSigOfNodes(myGOdataBP, score(resultFisherBP), firstSigNodes =n.nodes, useInfo = "all")
#showSigOfNodes(myGOdataBP, score(resultFisherBP), firstSigNodes =n.nodes, useInfo = "def")

setwd(paste0(path.to.save.files))
printGraph(myGOdataBP, resultFisherBP, 
           firstSigNodes = n.nodes, fn.prefix= paste0(feature.name,"_corr_BP_",up_down), useInfo="all", pdfSW= T)
#printGraph(myGOdataBP, resultFisherBP, 
#           firstSigNodes = n.nodes, fn.prefix= paste0(feature.name,"_corr_BP_pie_",up_down), useInfo="all", pdfSW= T, type = "type")
## ## ##
## ## ##
names(allResBP_view)[4] = "Significant Genes"
bubble_BP = allResBP_view %>% 
  ggplot(aes(Annotated, reorder(Term,Annotated), size = `Significant Genes`, color = weightFisher)) + 
  scale_color_gradient("weightFisher (<0.05)", low="#D14444", high="black") + #theme_classic() +
  labs(title = paste0("Biological Process - corr to ",feature.name," (",up_down,")"), x = "Annotated Genes", y = "") +
  theme(plot.title=element_text(hjust=0.5),
        legend.key.size = unit(0.45, 'cm'),
        legend.title = element_text(size=9.5)) + 
  geom_point()
ggsave(paste0(path.to.save.files,"/",feature.name,"_corr_BP_plot_",up_down,"_topGO.pdf"), 
       plot = bubble_BP, width = 5.24, height = 4.86)
## ## ##
## ## ##


#### CC Ontology ####
myGOdataCC <- new("topGOdata", ontology= "CC", allGenes=geneList, 
                  annot= annFUN.gene2GO, gene2GO= geneID2GO, nodeSize = 5)
sgCC <- sigGenes(myGOdataCC)
str(sgCC)
numSigGenes(myGOdataCC)

resultFisherCC <- runTest(myGOdataCC, algorithm = "weight01", statistic = "fisher")
allResCC <- GenTable(myGOdataCC, weightFisher= resultFisherCC, 
                     orderBy= "weightFisher", ranksOf= "weightFisher", topNodes= 200)
allResCC$weightFisher = as.numeric(allResCC$weightFisher)
allResCC_view = allResCC[allResCC$weightFisher <= 0.05,]

write.csv(allResCC_view, paste0(path.to.save.files,"/",feature.name,"_corr_CC_",up_down,"_topGO.csv"), quote = F, row.names = F)

#showSigOfNodes(myGOdataCC, score(resultFisherCC), firstSigNodes =n.nodes, useInfo = "all")
#showSigOfNodes(myGOdataCC, score(resultFisherCC), firstSigNodes =n.nodes, useInfo = "def")

printGraph(myGOdataCC, resultFisherCC, firstSigNodes = n.nodes, fn.prefix= paste0(feature.name,"_corr_CC_",up_down), 
           useInfo="all", pdfSW= T)
## ## ##
## ## ##
names(allResCC_view)[4] = "Significant Genes"
bubble_CC = allResCC_view %>% 
  ggplot(aes(Annotated, reorder(Term,Annotated), size = `Significant Genes`, color = weightFisher)) + 
  scale_color_gradient("weightFisher (<0.05)", low="#DFA377", high="black") + #theme_classic() +
  labs(title = paste0("Cellular Component - corr to ",feature.name," (",up_down,")"), x = "Annotated Genes", y = "") +
  theme(plot.title=element_text(hjust=0.5),
        legend.key.size = unit(0.45, 'cm'),
        legend.title = element_text(size=9.5)) + 
  geom_point()
ggsave(paste0(path.to.save.files,"/",feature.name,"_corr_CC_plot_",up_down,"_topGO.pdf"), 
       plot = bubble_CC, width = 5.24, height = 4.86)
## ## ##
## ## ##


### MF Ontolopy  ###
myGOdataMF <- new("topGOdata", ontology= "MF", allGenes=geneList, annot= annFUN.gene2GO,
                  gene2GO= geneID2GO, nodeSize= 5)
sgMF <- sigGenes(myGOdataMF)
str(sgMF)
numSigGenes(myGOdataMF)

resultFisherMF <- runTest(myGOdataMF, algorithm = "weight01", statistic = "fisher")
allResMF <- GenTable(myGOdataMF, weightFisher= resultFisherMF, orderBy= "weightFisher",
                     ranksOf= "weightFisher", topNodes= 200)
allResMF$weightFisher = as.numeric(allResMF$weightFisher)
allResMF_view = allResMF[allResMF$weightFisher <= 0.05,]

write.csv(allResMF_view, paste0(path.to.save.files,"/",feature.name,"_corr_MF_",up_down,"_topGO.csv"), quote = F, row.names = F)

#showSigOfNodes(myGOdataMF, score(resultFisherMF), firstSigNodes =n.nodes, useInfo = "all")
#showSigOfNodes(myGOdataMF, score(resultFisherMF), firstSigNodes =n.nodes, useInfo = "def")

printGraph(myGOdataMF, resultFisherMF, firstSigNodes = n.nodes, fn.prefix= paste0(feature.name,"_corr_MF_",up_down),
           useInfo="all", pdfSW= T)
#printGraph(myGOdataMF, resultFisherMF, firstSigNodes = n.nodes, fn.prefix= paste0(treatment,"_MF"),
#           useInfo="def", pdfSW= T)
## ## ##
## ## ##
names(allResMF_view)[4] = "Significant Genes"
bubble_MF = allResMF_view %>% 
  ggplot(aes(Annotated, reorder(Term,Annotated), size = `Significant Genes`, color = weightFisher)) + 
  scale_color_gradient("weightFisher (<0.05)", low="#2E4AA4", high="black") + #theme_classic() +
  labs(title = paste0("Molecular Function - corr to ",feature.name," (",up_down,")"), x = "Annotated Genes", y = "") +
  theme(plot.title=element_text(hjust=0.5),
        legend.key.size = unit(0.45, 'cm'),
        legend.title = element_text(size=9.5)) + 
  geom_point()
ggsave(paste0(path.to.save.files,"/",feature.name,"_corr_MF_plot_",up_down,"_topGO.pdf"), 
       plot = bubble_MF, width = 5.24, height = 4.86)
## ## ##
## ## ##


setwd(ori_path)
}
