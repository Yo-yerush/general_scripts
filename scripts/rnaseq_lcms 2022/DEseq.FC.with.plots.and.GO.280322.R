deseq_fc =  function(A.B_VS_c, experiment, samples.number, exp.treatment=NA, control=NA, group.A = NA,
                     with.GO = T, PCA.treatment.labels=F, PCA.id.labels=F, PCA.only=F, heatmap.plots=T, save.csv=T, save.plots=T, path=getwd()) {

  setwd(path)
  if (PCA.only == T) {
    save.csv = F
    save.plots = F
  } else {
    if (file.exists(paste(path,experiment,"/statistics/",A.B_VS_c, sep="/")) == T) {
      stop(paste("change <",A.B_VS_c,">, file is already exists", sep = ""))
    }
  }
  
  
  if (is.na(group.A)[1] == T) {
    
    if (is.na(exp.treatment) == T) {
      stop("experiment treatment name from coldata file: <exp.treatment> VS <control> (mut1 VS WT)")
    }
    if (is.na(control) == T) {
      stop("control treatment name from coldata file: <exp.treatment> VS <control> (mut1 VS WT)")
    }
    contrast = c("exp", exp.treatment, control)
    
  } else {
    contrast = c("group.for.stat", "A", "B")
  }
  
  ##################################
  logFC.up.file = paste(A.B_VS_c,"FC.up.DE.csv", sep=".")
  logFC.down.file = paste(A.B_VS_c,"FC.down.DE.csv", sep=".")
  norm.file = paste("norm", A.B_VS_c,"DE.csv", sep=".")
  
  pre_main_title = gsub("\\."," & ", A.B_VS_c)
  main_title = gsub("_", "  ", pre_main_title)
  
  #################################
  library(dplyr)
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(tximport)
  library("RColorBrewer")
  library(pheatmap)
  
  # coldata file
  setwd(paste(path,experiment ,sep="/"))
  info = read.table(paste("coldata",experiment,"txt", sep = "."), header = T, sep = '\t')
  info$sample = as.character(info$sample)
  
  if (is.na(group.A)[1] == T) {
    info = info[samples.number,]
  } else {
    info$group.for.stat = "B"
    info$group.for.stat[group.A] = "A"
    info = info[samples.number,]
  }
  
  files = c(paste(info$x,".genes.results", sep = ""))
  names(files) = c(info$x)
  
  
  ########################## DESeq2 ######################
  txData = tximport(files, 'rsem')
  txData$length[txData$length == 0] = 1
  if (contrast[1] == "exp") {
    dds = DESeqDataSetFromTximport(txData, info, ~exp)
  } else {
    dds = DESeqDataSetFromTximport(txData, info, ~group.for.stat)
  }
  
  keep = rowSums(counts(dds)) >=10 # remove low expressed genes
  dds = dds[keep,]
  
  if (PCA.only == F) {
    
    ddsDE = DESeq(dds)
    res = results(ddsDE, contrast = contrast, alpha = 0.05)
    
    samples.deseq = data.frame(rownames(res), res$log2FoldChange, res$padj)
    names(samples.deseq) = c("transcript_id", "log2FoldChange", "padj")
    
    summary(res)
    
    
    ################ merge names and samples by ACCESSION column, save csv file #############
    setwd(path)
    
    if (with.GO == T) {
      transcript_id.and.description.with.GO = read.csv("transcript_id.and.description.with.GO.csv")
      df_merge.pre = merge.data.frame(samples.deseq, transcript_id.and.description.with.GO, by = "transcript_id", all.x = T)
      df_merge.pre_2 = df_merge.pre[,c(1,4,2,3,5:8)]
    } else {
      # transcript_id and description from gtf file
      transcript_id.and.description = read.csv("transcript_id.and.description.csv", sep = '\t', header = F)
      names(transcript_id.and.description) = c("transcript_id", "description")
      
      df_merge.pre = merge.data.frame(samples.deseq, transcript_id.and.description, by = "transcript_id", all.x = T)
      df_merge.pre_2 = df_merge.pre[,c(1,4,2,3)]
    }
    
    df_merge.pre_3 = df_merge.pre_2[df_merge.pre_2$padj <= 0.05,] ## only p<0.05
    
    df_merge.pre_logFC.up = df_merge.pre_3[df_merge.pre_3$log2FoldChange >= 1,] ## log2FC > 1
    df_merge.pre_logFC.down = df_merge.pre_3[df_merge.pre_3$log2FoldChange <= -1,] ## log2FC < -1
    
    no.NA_up = df_merge.pre_logFC.up[complete.cases(df_merge.pre_logFC.up$padj), ]  ## remove NA
    no.NA_down = df_merge.pre_logFC.down[complete.cases(df_merge.pre_logFC.down$padj), ] ## remove NA
    
    df_merge_orderd_by_padj_up = no.NA_up[order(no.NA_up$padj, decreasing = T),]  ## orderd by padj
    df_merge_orderd_by_fc_up = no.NA_up[order(no.NA_up$log2FoldChange, decreasing = T),]  ## orderd by FC
    
    df_merge_orderd_by_padj_down = no.NA_down[order(no.NA_down$padj),]  ## orderd by padj
    df_merge_orderd_by_fc_down = no.NA_down[order(no.NA_down$log2FoldChange),]  ## orderd by FC
    
    #norm counts file
    normCounts = counts(ddsDE, normalized = T)
    
    if (save.csv == T) {
      dir.create(paste(path,experiment,"statistics",A.B_VS_c,sep="/"))
      setwd(paste(path,experiment,"statistics",A.B_VS_c,sep="/"))
      write.csv(df_merge_orderd_by_fc_up, logFC.up.file, row.names=FALSE)
      write.csv(df_merge_orderd_by_fc_down, logFC.down.file, row.names=FALSE)
      write.csv(df_merge.pre_2, paste(experiment,"all",A.B_VS_c,"DE.csv", sep="."), row.names=FALSE)
      write.csv(normCounts, norm.file)
    }
  }
  
  ###################### plots ###############################
  if (PCA.treatment.labels == T) {
    if (contrast[1] == "exp") {
      Label = info$exp
    } else {
      Label = info$group.for.stat
    } 
  } else {
    if (PCA.id.labels == T) {
      Label = info$x
    } else {
      Label = NA
    }
  }
  
  # plotPCA with ggplot
  vsd <- varianceStabilizingTransformation(dds)
  dds.est <- estimateSizeFactors(dds) # also possible to perform custom transformation
  se <- SummarizedExperiment(log2(counts(dds.est, normalized=TRUE) + 1), colData=colData(dds)) # shifted log of normalized counts
  data = plotPCA(DESeqTransform(se), intgroup = "sample", returnData = T)
  percentVar <- round(100 * attr(data, "percentVar"))
  plotPCA = ggplot(data, aes(PC1, PC2, color=sample)) +
    geom_point(size=3) +
    labs(title = main_title) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme_classic() + geom_hline(yintercept=0) + geom_vline(xintercept = 0) +
    geom_text(aes(label=Label),hjust=0.25, vjust=-0.5, show_guide = F) 
  
  if (PCA.only == F) {
    if (heatmap.plots ==T) {
      # heatmap
      file1 = data.frame(df_merge.pre_2$transcript_id)
      file1.2 = data.frame(df_merge.pre_3$transcript_id) # p < 0.05
      names(file1) = c("transcript_id")
      names(file1.2) = c("transcript_id")
      
      file2.0 = normCounts
      file2 = cbind.data.frame(row.names(file2.0), file2.0)
      file2$transcript_id = file2$`row.names(file2.0)`
      file2 = file2[,-1]
      
      merge1 = merge.data.frame(file1, file2, by = "transcript_id")
      merge1.2 = merge1[,-1]
      merge2 = merge.data.frame(file1.2, file2, by = "transcript_id")
      merge2.2 = merge2[,-1]
      
      plot.heatmap.F = pheatmap(merge1.2, scale = "row", show_rownames = F, main = main_title, cluster_cols = F)
      plot.heatmap.T = pheatmap(merge1.2, scale = "row", show_rownames = F, main = paste(main_title,",  clusterd", sep = ""), cluster_cols = T)
      plot.heatmap.sig.T = pheatmap(merge2.2, scale = "row", show_rownames = F, main = paste(main_title,",  clusterd by statistically significant (P=<0.05)", sep = ""), cluster_cols = T)
    }
    
    # volcano plot
    plot_volcano <- function(res_obj,
                             FDR = 0.05,
                             xlim_both_sides = NULL,
                             ylim_up = NULL,
                             vlines = NULL,
                             title = NULL,
                             intgenes = NULL,
                             intgenes_color = "steelblue",
                             labels_intgenes = TRUE,
                             labels_repel = TRUE) {
      mydf <- as.data.frame(res_obj)
      mydf$id <- rownames(mydf)
      mydf$isDE <- ifelse(is.na(res_obj$padj), FALSE, res_obj$padj < FDR)
      
      mydf <- mydf[mydf$baseMean > 0, ]
      
      p <- ggplot(mydf, aes_string(x = "log2FoldChange", y = "-log10(pvalue)")) +
        geom_point(aes_string(color = "isDE"), alpha = 0.4)
      
      
      if (!is.null(ylim_up)) {
        p <- p + coord_cartesian(ylim = c(0, ylim_up), xlim = c(-xlim_both_sides, xlim_both_sides))
      } else {
        p <- p + coord_cartesian(ylim = c(0, 20), xlim = c(-xlim_both_sides, xlim_both_sides))
      }
      
      if (!is.null(title)) {
        p <- p + ggtitle(title)
      }
      
      p <- p + theme_bw() +
        scale_colour_manual(
          name = paste0("FDR = ", FDR),
          values = c("black", "red"),
          labels = c("nonDE", "DE")
        )
      
      p <- p + geom_vline(aes(xintercept = 1), col = "lightblue", alpha = 0.6, size = 1.5) +
        geom_vline(aes(xintercept = -1), col = "lightblue", alpha = 0.6, size = 1.5)
      
      if (!is.null(intgenes)) {
        if ("symbol" %in% colnames(mydf)) {
          # use the gene names
          df_intgenes <- mydf[mydf$symbol %in% intgenes, ]
          df_intgenes$myids <- df_intgenes$symbol
        } else {
          # use whatever is there as id
          df_intgenes <- mydf[rownames(mydf) %in% intgenes, ]
          df_intgenes$myids <- rownames(df_intgenes)
        }
        
        # df_intgenes <- mydf[mydf$symbol %in% intgenes,]
        
        p <- p + geom_point(data = df_intgenes, aes_string("log2FoldChange", "-log10(pvalue)"), color = intgenes_color, size = 4)
        
        if (labels_intgenes) {
          
          if (labels_repel) {
            p <- p + geom_text_repel(
              data = df_intgenes, aes_string("log2FoldChange", "-log10(pvalue)", label = "myids"),
              color = intgenes_color, size = 5
            )
          } else {
            p <- p + geom_text(
              data = df_intgenes, aes_string("log2FoldChange", "-log10(pvalue)", label = "myids"),
              color = intgenes_color, size = 5, hjust = 0.25, vjust = -0.75
            )
          }
          
        }
      }
      
      p
    }
    
    
    
    setwd(path)
    vol.plot = plot_volcano(res, xlim_both_sides = 15, ylim_up = 30)
    vol.plot
    
    if (save.plots == T) {
      dir.create(paste(path,experiment,"statistics",A.B_VS_c,"plots",sep="/"))
      setwd(paste(path,experiment,"statistics",A.B_VS_c,"plots",sep="/"))
      # MA plot - 
      resApeT <- lfcShrink(ddsDE, coef=2, type="apeglm", lfcThreshold=1)
      plotMA(resApeT, ylim=c(-6,6), cex=.8) %>% abline(h=c(-1,1), col="dodgerblue", lwd=2)
      dev.copy(png, paste("MA", A.B_VS_c, "png", sep = "."))
      dev.off()
      ggsave(paste("PCA", A.B_VS_c, "png", sep = "."), plot = plotPCA)
      write.csv(data, paste("PCA", A.B_VS_c, "table","csv", sep = "."),row.names = F)
      ggsave(paste("heatmap", A.B_VS_c, "png", sep = "."), plot = plot.heatmap.F)
      ggsave(paste("heatmap.clusterd", A.B_VS_c, "png", sep = "."), plot = plot.heatmap.T)
      ggsave(paste("heatmap.sig.clusterd", A.B_VS_c, "png", sep = "."), plot = plot.heatmap.sig.T)
      ggsave(paste("volcano", A.B_VS_c, "png", sep = "."), plot = vol.plot)
    } else {
      resApeT <- lfcShrink(ddsDE, coef=2, type="apeglm", lfcThreshold=1)
      plotMA(resApeT, ylim=c(-6,6), cex=.8) %>% abline(h=c(-1,1), col="dodgerblue", lwd=2)
    }
  }
  plotPCA
}