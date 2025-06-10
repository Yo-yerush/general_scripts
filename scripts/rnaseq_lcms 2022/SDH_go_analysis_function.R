top.GO.fun = function(treatment = "high_VS_low",
log2FC_is_positive = c(T, 0.5), # true/false for positive/negetive and 0.5 got log2FC trashold
GO_Ontology_type, # BP, CC, MF
n.nodes = 10,
get_GO_term_genes = c(F, "GO_term"),
save.files = T) {


library(topGO)
library(Rgraphviz)
library(dplyr)

##########################################
  exp_1_file = read.csv("P:/yonatan/RNAseq_yonatan_2021/DESeq2/peel/statistics/all_high_VS_low/peel.all.all_high_VS_low.DE.csv")
##########################################
    
if (log2FC_is_positive[1] == T) {
  up_down = "up"
} else if (log2FC_is_positive[1] == F) {
  up_down = "down"
} 
  
start_path = getwd()

uniprot_pome = read.csv("P:/pomegranate RNAseq 2022/GO files/uniprot-pomegranate.csv")[,c(1,3)]
uniprot_2_pro = read.csv("P:/pomegranate RNAseq 2022/GO files/uni.by.XP.csv")[,-c(2:4)]
pro_2_GO = merge.data.frame(uniprot_2_pro, uniprot_pome, all.x = T)[,-1]
pro_2_GO = pro_2_GO[!pro_2_GO$Gene.ontology.IDs == "",]
pro_2_GO$Gene.ontology.IDs = gsub("; ",",",pro_2_GO$Gene.ontology.IDs)
pro_2_GO$protein_id = noquote(pro_2_GO$protein_id)

trnscript_GO_id = trnscript_GO_id_row
trnscript_GO_id = trnscript_GO_id[trnscript_GO_id$padj <= 0.05,]

if (log2FC_is_positive[1] == T) {
  trnscript_GO_id = trnscript_GO_id[trnscript_GO_id$log2FoldChange >= log2FC_is_positive[2],]
} else if (log2FC_is_positive[1] == F) {
  trnscript_GO_id = trnscript_GO_id[trnscript_GO_id$log2FoldChange <= -log2FC_is_positive[2],]
}

#trnscript_GO_id = na.omit(trnscript_GO_id)

diff_exp_pro = trnscript_GO_id[,2]

#####
merged.test = merge.data.frame(trnscript_GO_id_row, pro_2_GO, by = "protein_id")
merged.test = merged.test[,names(pro_2_GO)]
#####

#write.table(pro_2_GO, "P:/pomegranate RNAseq 2022/GO files/remove_1.txt", 
#            row.names = F, col.names = F,quote = 2, sep = "\t")
write.table(merged.test, 
            paste(path_for_GO_analysis_results_folder, "remove_1.txt", sep = ""), 
            row.names = F, col.names = F,quote = 2, sep = "\t")
write.table(diff_exp_pro,
            paste(path_for_GO_analysis_results_folder, "remove_2.txt", sep = ""),
            row.names = F, col.names = F, quote = F)

geneID2GO <- readMappings(file = paste(path_for_GO_analysis_results_folder, "remove_1.txt",
                                       sep = ""))   
geneUniverse <- names(geneID2GO)   
genesOfInteres <- read.table(file = paste(path_for_GO_analysis_results_folder, "remove_2.txt", 
                                          sep = ""), header = F)
genesOfInteres <- as.character(genesOfInteres$V1)
geneList<- factor(as.integer(geneUniverse %in% genesOfInteres))
names(geneList) <- geneUniverse

file.remove(c(paste(path_for_GO_analysis_results_folder, "remove_1.txt", sep = ""),
            paste(path_for_GO_analysis_results_folder, "remove_2.txt", sep = "")))



############################################################
if (GO_Ontology_type == "BP") {
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

if (save.files == T) {
  if (dir.exists(paste(path_for_GO_analysis_results_folder,"_", treatment,
                       sep = "")) == F) {
      dir.create(paste(path_for_GO_analysis_results_folder,"_", treatment, sep = ""))
  }
  write.csv(allResBP_view, paste(path_for_GO_analysis_results_folder,"_", treatment,
                            "/","_", treatment,  ".",
                            GO_Ontology_type ,"." ,up_down,".topGO.csv", sep = ""),
            quote = F, row.names = F)
}

showSigOfNodes(myGOdataBP, score(resultFisherBP), firstSigNodes =n.nodes, useInfo = "all")
showSigOfNodes(myGOdataBP, score(resultFisherBP), firstSigNodes =n.nodes, useInfo = "def")

setwd(paste(path_for_GO_analysis_results_folder,"_", treatment, sep = ""))
printGraph(myGOdataBP, resultFisherBP, 
           firstSigNodes = n.nodes, fn.prefix= paste0("_",treatment,"_BP_",up_down), useInfo="all", pdfSW= T)
#printGraph(myGOdataBP, resultFisherBP, 
#           firstSigNodes = 5, fn.prefix= paste0(treatment,"_BP"), useInfo="def", pdfSW= T)



####################################################
} else if (GO_Ontology_type == "CC") {
#### CC Ontology ####
myGOdataCC <- new("topGOdata", ontology= "CC", allGenes=geneList, 
                  annot= annFUN.gene2GO, gene2GO= geneID2GO, nodeSize = 5)
sgCC <- sigGenes(myGOdataCC)
str(sgCC)
numSigGenes(myGOdataCC)


resultFisherCC <- runTest(myGOdataCC, algorithm = "weight01", statistic = "fisher")
allResCC <- GenTable(myGOdataCC, weightFisher= resultFisherCC, 
                     orderBy= "weightFisher", ranksOf= "weightFisher", topNodes=length(resultFisherCC@score))
allResCC$weightFisher = as.numeric(allResCC$weightFisher)

allResCC_view = allResCC[allResCC$weightFisher <= 0.05,]

if (save.files == T) {
  if (dir.exists(paste(path_for_GO_analysis_results_folder,"_", treatment,
                       sep = "")) == F) {
      dir.create(paste(path_for_GO_analysis_results_folder,"_", treatment, sep = ""))
  }
  write.csv(allResCC_view, paste(path_for_GO_analysis_results_folder,"_", treatment,
                            "/","_", treatment,  ".",
                            GO_Ontology_type ,"." ,up_down,".topGO.csv", sep = ""),
            quote = F, row.names = F)
}

showSigOfNodes(myGOdataCC, score(resultFisherCC), firstSigNodes =n.nodes, useInfo = "all")
showSigOfNodes(myGOdataCC, score(resultFisherCC), firstSigNodes =n.nodes, useInfo = "def")

setwd(paste(path_for_GO_analysis_results_folder,"_", treatment, sep = ""))
printGraph(myGOdataCC, resultFisherCC, firstSigNodes = n.nodes, fn.prefix= paste0("_",treatment,"_CC_",up_down), 
           useInfo="all", pdfSW= T)
#printGraph(myGOdataCC, resultFisherCC, firstSigNodes = 5, fn.prefix= paste0(treatment,"_CC"), 
#           useInfo="def", pdfSW= T)




#########################################################
} else if (GO_Ontology_type == "MF") {
### MF Ontolopy  ###
myGOdataMF <- new("topGOdata", ontology= "MF", allGenes=geneList, annot= annFUN.gene2GO,
                  gene2GO= geneID2GO, nodeSize= 5)
sgMF <- sigGenes(myGOdataMF)
str(sgMF)
numSigGenes(myGOdataMF)


resultFisherMF <- runTest(myGOdataMF, algorithm = "weight01", statistic = "fisher")
allResMF <- GenTable(myGOdataMF, weightFisher= resultFisherMF, orderBy= "weightFisher",
                     ranksOf= "weightFisher", topNodes=length(resultFisherMF@score))
allResMF$weightFisher = as.numeric(allResMF$weightFisher)

allResMF_view = allResMF[allResMF$weightFisher <= 0.05,]

if (save.files == T) {
  if (dir.exists(paste(path_for_GO_analysis_results_folder,"_", treatment,
                       sep = "")) == F) {
      dir.create(paste(path_for_GO_analysis_results_folder,"_", treatment, sep = ""))
  }
  write.csv(allResMF_view, paste(path_for_GO_analysis_results_folder,"_", treatment,
                            "/","_", treatment,  ".",
                            GO_Ontology_type ,"." ,up_down,".topGO.csv", sep = ""),
            quote = F, row.names = F)
}

showSigOfNodes(myGOdataMF, score(resultFisherMF), firstSigNodes =n.nodes, useInfo = "all")
showSigOfNodes(myGOdataMF, score(resultFisherMF), firstSigNodes =n.nodes, useInfo = "def")

setwd(paste(path_for_GO_analysis_results_folder,"_", treatment, sep = ""))
printGraph(myGOdataMF, resultFisherMF, firstSigNodes = n.nodes, fn.prefix= paste0("_",treatment,"_MF_",up_down),
           useInfo="all", pdfSW= T)
#printGraph(myGOdataMF, resultFisherMF, firstSigNodes = 5, fn.prefix= paste0(treatment,"_MF"),
#           useInfo="def", pdfSW= T)
}
 
if (get_GO_term_genes[1] == T) {
  # Extract all sig genes from GO 
  myterms <- get_GO_term_genes[2]
  mygenes <- genesInTerm(myGOdataBP, myterms)
  
  for (i in 1:length(myterms))
  {
    myterm <- myterms[i]
    mygenesforterm <- mygenes[myterm][[1]]
    mygenesforterm <- paste(mygenesforterm, collapse=',')
    print(paste("Term: ",myterm," genes:",mygenesforterm)) 
  }
  
  for (GoTerm in myterms) {
    print(GoTerm)
    a = unlist(mygenes[[GoTerm]])
    print(a)
    a_df = data.frame(gene_id = as.vector(a), 
                      GO_id = replicate(length(as.vector(a)), GoTerm))
    names(a_df)[1] = "protein_id"
  }
  
  tmp.t = merge.data.frame(trnscript_GO_id, a_df, by = "protein_id")
  tmp.t = tmp.t[,c(1,9,2:8)]
  
  if (GO_Ontology_type == "BP") {
      GO_Ontology_type_save = allResBP
  } else   if (GO_Ontology_type == "CC") {
    GO_Ontology_type_save = allResCC
  } else   if (GO_Ontology_type == "MF") {
    GO_Ontology_type_save = allResMF
  }
  
  write.csv(tmp.t, paste(path_for_GO_analysis_results_folder,"_", treatment,
                            "/","_", treatment,  ".sig.genes.",
                            GO_Ontology_type, ".", up_down, ".",
                         GO_Ontology_type_save[GO_Ontology_type_save$GO.ID == get_GO_term_genes[2],2],".csv", sep = ""),
            quote = F, row.names = F)
}
setwd(start_path)
}
