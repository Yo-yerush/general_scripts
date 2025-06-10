top.GO.fun = function(treatment,
                      log2FC_is_positive = c(T, 0.5), # true/false for positive/negetive and 0.5 got log2FC trashold
                      experiment,
                      GO_Ontology_type, # BP, CC, MF
                      pvalue = 0.05,
                      GO_algorithm = "classic", # classic, weight, weight01
                      path_for_GO_analysis_results_folder,
                      n.nodes = NULL,
                      get_GO_term_genes = c(F, "GO_term"),
                      save.files = T,
                      dotPlot = F,
                      width = 5,
                      height = 6) {
  
  
  library(topGO)
  library(Rgraphviz)
  library(dplyr)
  library(ggplot2)
  
  ##########################################
  exp_1_file = read.csv("P:/yonatan/RNAseq_yonatan_2021/DESeq2/peel/statistics/all_high_VS_low/peel.all.all_high_VS_low.DE.csv")
  exp_2_file = read.csv("P:/yonatan/RNAseq_yonatan_2021/DESeq2/clone/statistics/BanG.101_VS_DanA.116/clone.all.BanG.101_VS_DanA.116.DE.csv")
  exp_3_file = read.csv("P:/yonatan/RNAseq_yonatan_2021/DESeq2/transgenic/statistics/bHLH94-like_VS_EV/transgenic.all.bHLH94-like_VS_EV.DE.csv")
  
  exp_1_name = "peel"
  exp_2_name = "clone"
  exp_3_name = "transgenic"
  ##########################################
  
  if (log2FC_is_positive[1] == T) {
    up_down = "up"
  } else if (log2FC_is_positive[1] == F) {
    up_down = "down"
  } 
  
  start_path = getwd()
  
  uniprot_pome = read.csv("P:/pomegranate/pomegranate RNAseq 2022/GO files/uniprot-pomegranate.csv")[,c(1,3)]
  uniprot_2_pro = read.csv("P:/pomegranate/pomegranate RNAseq 2022/GO files/uni.by.XP.csv")[,-c(2:4)]
  pro_2_GO = merge.data.frame(uniprot_2_pro, uniprot_pome, all.x = T)[,-1]
  pro_2_GO = pro_2_GO[!pro_2_GO$Gene.ontology.IDs == "",]
  pro_2_GO$Gene.ontology.IDs = gsub("; ",",",pro_2_GO$Gene.ontology.IDs)
  pro_2_GO$protein_id = noquote(pro_2_GO$protein_id)
  
  if (experiment == exp_1_name) {
    trnscript_GO_id_row = exp_1_file
  } else if (experiment == exp_2_name) {
    trnscript_GO_id_row = exp_2_file
  } else if (experiment == exp_3_name) {
    trnscript_GO_id_row = exp_3_file
  }
  
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
  #### type Ontology  ####
  myGOdata <- new("topGOdata", ontology= GO_Ontology_type, allGenes=geneList, 
                  annot= annFUN.gene2GO, gene2GO= geneID2GO, nodeSize = 5)
  
  sg <- sigGenes(myGOdata)
  str(sg)
  numSigGenes(myGOdata)

  resultFisher <- runTest(myGOdata, algorithm = GO_algorithm, statistic = "fisher")
  
  #resultclassic <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")
  #resultWeight <- runTest(myGOdata, algorithm = "weight", statistic = "fisher")
  #resultWeight01 <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
  
  
  allRes <- GenTable(myGOdata, Fisher = resultFisher,
                     orderBy= "Fisher", ranksOf= "Fisher", topNodes=length(resultFisher@score))
  allRes$Fisher = as.numeric(allRes$Fisher)
  allRes$Term = gsub(",", ";", allRes$Term)
  
  allRes_view = allRes[allRes$Fisher <= pvalue,]
  
  if (save.files == T) {
    if (dir.exists(paste0(path_for_GO_analysis_results_folder,experiment,"_", treatment)) == F) {
      dir.create(paste0(path_for_GO_analysis_results_folder,experiment,"_", treatment))
    }
    write.csv(allRes_view, paste0(path_for_GO_analysis_results_folder,experiment,"_", treatment,
                                  "/",experiment,"_", treatment,  ".",
                                  GO_Ontology_type ,"." ,up_down,".topGO.csv"),
              quote = F, row.names = F)
  }
  
  
  if (is.null(n.nodes)) {
    if (nrow(allRes_view) > 10) {
      n.nodes = 10
    } else {
      n.nodes = nrow(allRes_view)
    }
  }
  #showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes =n.nodes, useInfo = "all")
  #showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes =n.nodes, useInfo = "def")
  setwd(paste0(path_for_GO_analysis_results_folder,experiment,"_", treatment))
  printGraph(myGOdata, resultFisher, 
             firstSigNodes = n.nodes, fn.prefix= paste(experiment,treatment,GO_Ontology_type,up_down,sep = "_"),
             useInfo="all", pdfSW= T)
  #printGraph(myGOdata, resultFisher, 
  #           firstSigNodes = 5, fn.prefix= paste0(treatment,"_"), useInfo="def", pdfSW= T)
  
  
  
  ####################################################
  
  
  if (get_GO_term_genes[1] == T) {
    # Extract all sig genes from GO 
    myterms <- get_GO_term_genes[2]
    mygenes <- genesInTerm(myGOdata, myterms)
    
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
    
    
    GO_Ontology_type_save = allRes
    
    write.csv(tmp.t, paste(path_for_GO_analysis_results_folder,experiment,"_", treatment,
                           "/",experiment,"_", treatment,  ".sig.genes.",
                           GO_Ontology_type, ".", up_down, ".",
                           GO_Ontology_type_save[GO_Ontology_type_save$GO.ID == get_GO_term_genes[2],2],".csv", sep = ""),
              quote = F, row.names = F)
  }
  
  if (dotPlot) {
    if (GO_Ontology_type == "BP") {
      low_col = "#cf534c"
      type_name = "Biological Process"
    } else if (GO_Ontology_type == "CC") {
      low_col = "#d1973f"
      type_name = "Cellular Component"
    } else if (GO_Ontology_type == "MF") {
      low_col = "#6397eb"
      type_name = "Molecular Function"
    }
    #names(allRes_view)[4] = "Significant Genes"
    
    bubble = allRes_view %>% 
      ggplot(aes(Significant, reorder(Term,Significant), size = Annotated, color = Fisher)) + 
      scale_color_gradient("p.value", low=low_col, high="black") + #theme_classic() +
      labs(#title = paste0(type_name," - ",up_down,"regulated transcripts"),
           x = "Significant", y = "") + theme_bw() + 
      theme(
        #plot.title=element_text(hjust=0.5),
            legend.key.size = unit(0.25, 'cm'),
            legend.title = element_text(size=9.5),
            text = element_text(family = "serif")) + 
      geom_point()
    ggsave(paste0(path_for_GO_analysis_results_folder,experiment,"_", treatment,
                  "/",treatment,"_",GO_Ontology_type,"_dotPlot_",up_down,"_topGO.png"), 
           plot = bubble, width = width, height = height)
  }
  
  setwd(start_path)
}
