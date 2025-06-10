library(DMRcaller)
library(GenomicFeatures)
library(dplyr)
library(parallel)

###################################
# upload CX files
var1_path <- c(
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S18/methylation_extractor/S18_R1_bismark_bt2_pe.CX_report.txt",
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S19/methylation_extractor/S19_R1_bismark_bt2_pe.CX_report.txt"
)

var2_path <- c(
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S9/methylation_extractor/S9_R1_bismark_bt2_pe.CX_report.txt",
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S10/methylation_extractor/S10_R1_bismark_bt2_pe.CX_report.txt",
  "/home/yoyerush/yo/methylome_pipeline/Bismark/res_310523/S11/methylation_extractor/S11_R1_bismark_bt2_pe.CX_report.txt"
)

var1 <- mclapply(var1_path, readBismark, mc.cores = 2)
var1_pool <- poolMethylationDatasets(GRangesList(var1))

var2 <- mclapply(var2_path, readBismark, mc.cores = 3)
var2_pool <- poolMethylationDatasets(GRangesList(var2))

var1_pool <- rename_seq(var1_pool)
var2_pool <- rename_seq(var2_pool)

###################################

ann_file <- read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/annotation_files/Methylome.At_annotations.csv.gz") %>%
  select(-width)

TE_file = read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/annotation_files/TAIR10_Transposable_Elements.txt", sep = "\t") %>%
  mutate(seqnames = NA) %>% # Add a new column with NA values
  mutate(type = "transposable_element", gene_model_type = "transposable_element") %>%
  dplyr::select(seqnames, Transposon_min_Start, Transposon_max_End, orientation_is_5prime, type, Transposon_Name, gene_model_type) %>%
  dplyr::rename(gene_id = Transposon_Name)
  for (i in 1:5) {
    TE_file$seqnames[grep(paste0("AT", i, "TE"), TE_file$gene_id)] = paste0("Chr", i)
  }
  TE_file$orientation_is_5prime = gsub("true", "+", TE_file$orientation_is_5prime)
  TE_file$orientation_is_5prime = gsub("false", "-", TE_file$orientation_is_5prime)
  names(TE_file)[1:4] = c("seqnames", "start", "end", "strand")


gff3_file = rbind(ann_file, TE_file) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

#gff3_file <- rtracklayer::import("/home/yoyerush/yo/TAIR10/TAIR10 gff3/TAIR10_GFF3_genes_transposons.gff")

###################################

genePlot_fun = function(tair_id, gff3 = gff3_file, var1_pool_f = var1_pool, var2_pool_f = var2_pool) {

# get gene position from the gff3 file
gene_gr <- as.data.frame(gff3) %>%
  filter(type == "gene") %>%
  filter(gene_id == tair_id) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

if (as.character(strand(gene_gr)) == "+") {
  start(gene_gr) <- start(gene_gr) - 4000
  end(gene_gr) <- end(gene_gr) + 4000
} else {
  start(gene_gr) <- start(gene_gr) - 4000
  end(gene_gr) <- end(gene_gr) + 4000
}

chr_name <- as.character(gene_gr@seqnames@values[1])
start_pos <- start(gene_gr)
end_pos <- end(gene_gr)


# Filter var_pool by the positions of chr_name, start_pos, and end_pos
filtered_var1_pool <- var1_pool_f[seqnames(var1_pool_f) == chr_name & start(var1_pool_f) >= start_pos & end(var1_pool_f) <= end_pos]
filtered_var2_pool <- var2_pool_f[seqnames(var2_pool_f) == chr_name & start(var2_pool_f) >= start_pos & end(var2_pool_f) <= end_pos]
###################################

###################################
# DMRs file
columnNames <- c("seqnames", "start", "end", "width", "strand", "sumReadsM1", "sumReadsN1", "proportion1", "sumReadsM2", "sumReadsN2", "proportion2", "cytosinesCount", "context", "direction", "pValue", "regionType", "gene_id", "Symbol")
# DMRs = rbind(
cg <- rbind(read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/mto1_vs_wt/genome_annotation/CG/Genes_CG_genom_annotations.csv"),
read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/mto1_vs_wt/genome_annotation/CG/Promoters_CG_genom_annotations.csv")
) %>%
  select(all_of(columnNames)) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

chg <- rbind(read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/mto1_vs_wt/genome_annotation/CHG/Genes_CHG_genom_annotations.csv"),
read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/mto1_vs_wt/genome_annotation/CHG/Promoters_CHG_genom_annotations.csv")
) %>%
  select(all_of(columnNames)) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

chh <- rbind(read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/mto1_vs_wt/genome_annotation/CHH/Genes_CHH_genom_annotations.csv"),
read.csv("/home/yoyerush/yo/methylome_pipeline/Methylome.At/results/mto1_vs_wt/genome_annotation/CHH/Promoters_CHH_genom_annotations.csv")
) %>%
  select(all_of(columnNames)) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
# ) %>%

# DMRsList <- list("DMRs" = DMRs)
DMRsList <- list("CG" = cg, "CHG" = chg, "CHH" = chh)

###################################

id2symbol = data.frame(
  id = c(cg$gene_id, chg$gene_id, chh$gene_id),
  symbol = c(cg$Symbol, chg$Symbol, chh$Symbol)
) %>%
  filter(id == tair_id) %>%
  distinct(., .keep_all = TRUE)

main_title = ifelse(
  id2symbol$symbol == "",
  id2symbol$id,
  paste0(id2symbol$id, " (", id2symbol$symbol, ")")
)
###################################

dir.create("/home/yoyerush/yo/genePlot_TEs_n_gene_201124", showWarnings = FALSE)

svg(paste0("/home/yoyerush/yo/genePlot_TEs_n_gene_201124/genePlot_TEnGene_", main_title, ".svg"), width = 4.75, height = 2.75, family = "serif")
plotLocalMethylationProfile_yo(filtered_var1_pool,
  filtered_var2_pool,
  gene_gr,
  DMRsList,
  conditionsNames = c("WT", "mto1"),
  gff3,
  windowSize = 1000,
  context = c("CG", "CHG", "CHH"),
  main = main_title,
  col = c(
    "gray25", "#bf6828", "#000000", "#afad43",
    "#009E73", "#0072B2", "#CC79A7", "#D55E00", "#999999"
  )
)
dev.off()
}

##########################################



##########################################

genePlot_fun(tair_id = "AT3G01120")
genePlot_fun(tair_id = "AT1G48360")
genePlot_fun(tair_id = "AT1G42680")
genePlot_fun(tair_id = "AT4G00360")
genePlot_fun(tair_id = "AT1G53480")
genePlot_fun(tair_id = "AT5G35490")
genePlot_fun(tair_id = "AT4G00970")
genePlot_fun(tair_id = "AT5G35630")
genePlot_fun(tair_id = "AT2G15890")


##########################################
##########################################
### functions
rename_seq <- function(gr_obj) {
  gr_obj <- renameSeqlevels(gr_obj, gsub("NC_003070.9", "Chr1", seqlevels(gr_obj)))
  gr_obj <- renameSeqlevels(gr_obj, gsub("NC_003071.7", "Chr2", seqlevels(gr_obj)))
  gr_obj <- renameSeqlevels(gr_obj, gsub("NC_003074.8", "Chr3", seqlevels(gr_obj)))
  gr_obj <- renameSeqlevels(gr_obj, gsub("NC_003075.7", "Chr4", seqlevels(gr_obj)))
  gr_obj <- renameSeqlevels(gr_obj, gsub("NC_003076.8", "Chr5", seqlevels(gr_obj)))
  return(gr_obj)
}

.joinMethylationData <- function(cx1, cx2){
  
  overlaps <- findOverlaps(cx1, cx2)
  indexes <- which(!duplicated(queryHits(overlaps)))
  methylData <- GRanges(seqnames = seqnames(cx1[queryHits(overlaps)[indexes]]), 
                        ranges   = ranges(cx1[queryHits(overlaps)[indexes]]), 
                        strand   = strand(cx1[queryHits(overlaps)[indexes]]), 
                        context  = cx1$context[queryHits(overlaps)[indexes]],
                        trinucleotide_context = cx1$trinucleotide_context[queryHits(overlaps)[indexes]], 
                        readsM1  = cx1$readsM[queryHits(overlaps)[indexes]], 
                        readsN1  = cx1$readsN[queryHits(overlaps)[indexes]],                                   
                        readsM2  = cx2$readsM[subjectHits(overlaps)[indexes]], 
                        readsN2  = cx2$readsN[subjectHits(overlaps)[indexes]])
  return(methylData)
}

.movingAverage <- function(minPos, maxPos, pos, val, weights=1, windowSizeHalf=150, normalize = FALSE, kernelFunction="triangular", lambda=0.5) {
  
  if(length(weights) < length(val)){
    weights <- rep(weights, length.out=length(val))
  } else if(length(weights) > length(val)){
    weights <- weights[1:length(val)]
  }
  
  # Filter out NAs
  keepIndexes <- which(!is.na(pos) & !is.na(val) & !is.na(weights))
  pos <- pos[keepIndexes]
  weights <- weights[keepIndexes]  
  val <- val[keepIndexes]
    
  #set the values
  rawVector <- rep(0, maxPos - minPos + 2*windowSizeHalf + 1)
  rawVector[(pos - minPos + windowSizeHalf + 1)] <-weights*val
  
  normVector <- rep(0, maxPos - minPos + 2*windowSizeHalf + 1)
  normVector[(pos - minPos + windowSizeHalf + 1)] <-weights
  
  # Define the (triangular) kernel.
  if(kernelFunction =="uniform"){
    kernel <- rep(1, times=(2*windowSizeHalf + 1))  
  } else if(kernelFunction =="triangular"){
    kernel <- 1 - abs(-windowSizeHalf:windowSizeHalf)/windowSizeHalf 
  } else if(kernelFunction =="gaussian"){
    kernel <- .gaussianKernel(lambda, - windowSizeHalf:windowSizeHalf)
  } else if(kernelFunction =="epanechnicov"){
    kernel <- .epanechnicovKernel(-windowSizeHalf:windowSizeHalf, windowSizeHalf)
  } else{
    stop(paste("Unknown kernel function: ", kernelFunction, ". 
               It should be one of \"uniform\", \"triangular\", \"gaussian\", \"epanechnicov\"",sep=""))
  }
  
  kernel <- kernel / sum(kernel)
  
  if(windowSizeHalf >= 1){
    smoothedVector <- RcppRoll::roll_sum(rawVector, length(kernel), weights = kernel, normalize = normalize) / RcppRoll::roll_sum(normVector, length(kernel), weights = kernel, normalize = normalize)
  } else{
    smoothedVector <- rawVector 
  }
  
  return(smoothedVector)
}


.plotGeneticElements <- function(gff, region, col){

  seqname <- seqnames(region)
  minPos <- start(region)
  maxPos <- end(region)


  # Select the genes that lie in the region of interest.
  gff <- gff[queryHits(findOverlaps(gff, region))]
  # Chop off the ends of anything sticking out...
  start(gff) <- pmax(start(gff), minPos)
  end(gff)   <- pmin(end(gff), maxPos)

  genes <- gff[gff$type == 'gene']
  genesPos <- genes[strand(genes) == '+' | strand(genes) == '*']
  genesNeg <- genes[strand(genes) == '-' | strand(genes) == '*']
  exons <- gff[gff$type == 'exon']
  exons <- exons[overlapsAny(exons, genes)];
  exonsPos <- exons[strand(exons) == '+' | strand(exons) == '*']
  exonsNeg <- exons[strand(exons) == '-' | strand(exons) == '*']

  transposons <- gff[gff$type == 'transposable_element'];
  transposonsPos <- transposons[strand(transposons) == '+' | strand(transposons) == '*'];
  transposonsNeg <- transposons[strand(transposons) == '-' | strand(transposons) == '*'];

  negativeStrandPosition <- -0.175
  positiveStrandPosition <- -0.075


  #text(maxPos + (maxPos-minPos)/100, positiveStrandPosition, 'Sense', cex = 0.5);
  #text(maxPos + (maxPos-minPos)/100, negativeStrandPosition, 'Anti', cex = 0.5);
  text(maxPos + (maxPos-minPos)/100, positiveStrandPosition, '+', font = 2);
  text(maxPos + (maxPos-minPos)/100, negativeStrandPosition, '-', font = 2);
  lines(c(minPos, maxPos),c(-0.14,-0.14), lty=1, lwd=0.75, col="black")

  if(length(genesPos)>0) {
    segments(start(genesPos), positiveStrandPosition, end(genesPos), positiveStrandPosition);
    text(start(genesPos), -0.115, genesPos$ID, pos=4, cex=0.5);
  }
  if(length(genesNeg)>0) {
    segments(start(genesNeg), -0.175, end(genesNeg), negativeStrandPosition);
    text(start(genesNeg), -0.23, genesNeg$ID, pos=4, cex=0.5);
  }
  if(length(exonsPos)>0) {
    rect(start(exonsPos), -0.05, end(exonsPos), -0.09, col=col[1], border = NA);
  }
  if(length(exonsNeg)>0) {
    rect(start(exonsNeg), -0.16, end(exonsNeg), -0.2, col=col[1], border = NA);
  }
  if(length(transposonsPos)>0) {
    rect(start(transposonsPos), -0.05, end(transposonsPos), -0.09, col=col[2], border = col[2], density=30, angle=30);
    text(start(transposonsPos), -0.115, transposonsPos$ID, pos=4, cex=0.5);
  }
  if(length(transposonsNeg)>0) {
    rect(start(transposonsNeg), -0.16, end(transposonsNeg), -0.2, col=col[2], border = col[2], density=30, angle=30);
    text(start(transposonsNeg), -0.23, transposonsNeg$ID, pos=4, cex=0.5);
  }

}

plotLocalMethylationProfile_yo = function (methylationData1, methylationData2, region, DMRs = NULL, 
    conditionsNames = NULL, gff = NULL, windowSize = 150, context = "CG", 
    labels = NULL, col = NULL, main = "", plotMeanLines = TRUE, 
    plotPoints = TRUE) 
{
   
    numberOfConditions <- 2
    if (is.null(conditionsNames) | length(conditionsNames) < 
        numberOfConditions | !all(is.character(conditionsNames))) {
        conditionsNames <- paste("condition ", (1:numberOfConditions), 
            sep = "")
    }
   
    numberOfDMRs <- 0
    if (!is.null(DMRs)) {
        numberOfDMRs <- length(DMRs)
    }
    if (!is.null(labels) & (length(labels) < 1 | !is.character(labels))) {
        labels <- LETTERS[1:length(labels)]
    }

        cond1Color <- col[1]
        cond2Color <- col[2]
        geneColor <- col[3]
        TEColor <- col[4]
        DMRsColor <- col[5:length(col)]

    seqname <- seqnames(region)
    minPos <- start(region)
    maxPos <- end(region)
    methylationData <- .joinMethylationData(methylationData1, 
        methylationData2)
    hits <- findOverlaps(methylationData, region)
    localMethylationData <- methylationData[queryHits(hits)]
    contextMethylationData <- localMethylationData[localMethylationData$context %in% 
        context]
    ramp1 <- colorRampPalette(c("white", cond1Color))
    colramp1 <- ramp1(100)
    ramp2 <- colorRampPalette(c("white", cond2Color))
    colramp2 <- ramp2(100)
    proportion1 <- rep(0, length(contextMethylationData))
    index <- which(contextMethylationData$readsM1 >= 0 & contextMethylationData$readsN1 > 
        0)
    proportion1[index] <- contextMethylationData$readsM1[index]/contextMethylationData$readsN1[index]
    proportion2 <- rep(0, length(contextMethylationData))
    index <- which(contextMethylationData$readsM2 >= 0 & contextMethylationData$readsN2 > 
        0)
    proportion2[index] <- contextMethylationData$readsM2[index]/contextMethylationData$readsN2[index]
    maxColor <- max(c(contextMethylationData$readsN1[!is.na(contextMethylationData$readsN1)], 
        contextMethylationData$readsN2[!is.na(contextMethylationData$readsN2)]))
    par(mar = c(4, 4, 0, 1) + 0.1)
    plot(start(contextMethylationData), proportion1, col = colramp1[round(99 * 
        log(contextMethylationData$readsN1)/log(maxColor)) + 
        1], pch = 16, cex = 0.6, xlim = c(minPos, maxPos), ylim = c(-0.2, 
        1.2 + numberOfDMRs * 0.1), xlab = "genomic coordinate", ylab = "methylation proportion", 
        yaxt = "n", main = NULL, type = "n")
    axis(2, c(0, 0.5, 1))
    if (plotPoints) {
        points(start(contextMethylationData), proportion1, col = colramp1[round(99 * 
            log(contextMethylationData$readsN1)/log(maxColor)) + 
            1], pch = 16, cex = 0.6)
        points(start(contextMethylationData), proportion2, col = colramp2[round(99 * 
            log(contextMethylationData$readsN2)/log(maxColor)) + 
            1], pch = 16, cex = 0.6)
        if (!plotMeanLines) {
            #legend("topright", bty = "n", col = c(cond1Color, 
            #    cond2Color), legend = conditionsNames, pch = c(16, 
            #    16), horiz = TRUE)
        }
    }
    windowSizeHalf <- floor((windowSize - 1)/2)
    if (plotMeanLines) {
        positions <- start(region):end(region)
        movingAverageMethylReads1 <- .movingAverage(start(region), 
            end(region), start(contextMethylationData), contextMethylationData$readsM1, 
            windowSizeHalf = windowSizeHalf)
        movingAverageTotalReads1 <- .movingAverage(start(region), 
            end(region), start(contextMethylationData), contextMethylationData$readsN1, 
            windowSizeHalf = windowSizeHalf)
        movingAverageProportion1 <- movingAverageMethylReads1/movingAverageTotalReads1
        movingAverageMethylReads2 <- .movingAverage(start(region), 
            end(region), start(contextMethylationData), contextMethylationData$readsM2, 
            windowSizeHalf = windowSizeHalf)
        movingAverageTotalReads2 <- .movingAverage(start(region), 
            end(region), start(contextMethylationData), contextMethylationData$readsN2, 
            windowSizeHalf = windowSizeHalf)
        movingAverageProportion2 <- movingAverageMethylReads2/movingAverageTotalReads2
        maxColor <- max(c(movingAverageTotalReads1[!is.na(movingAverageTotalReads1)], 
            movingAverageTotalReads2[!is.na(movingAverageTotalReads2)]))
        bufferIndex <- !is.na(movingAverageProportion1) & !is.na(movingAverageTotalReads1)
        if (any(bufferIndex)) {
            minID <- min(which(bufferIndex))
            maxID <- max(which(bufferIndex))
            colorId <- round(99 * log(round(movingAverageTotalReads1))/log(maxColor)) + 
                1
            colorId <- round((colorId[minID:(maxID - 1)] + colorId[(minID + 
                1):(maxID)])/2)
            colorId[is.na(colorId)] <- 1
            segments(positions[minID:(maxID - 1)], movingAverageProportion1[minID:(maxID - 
                1)], positions[(minID + 1):(maxID)], movingAverageProportion1[(minID + 
                1):(maxID)], col = colramp1[colorId], lty = 1, 
                lwd = 2)
        }
        bufferIndex <- !is.na(movingAverageProportion2) & !is.na(movingAverageTotalReads2)
        if (any(bufferIndex)) {
            minID <- min(which(bufferIndex))
            maxID <- max(which(bufferIndex))
            colorId <- round(99 * log(round(movingAverageTotalReads2))/log(maxColor)) + 
                1
            colorId <- round((colorId[minID:(maxID - 1)] + colorId[(minID + 
                1):(maxID)])/2)
            colorId[is.na(colorId)] <- 1
            segments(positions[minID:(maxID - 1)], movingAverageProportion2[minID:(maxID - 
                1)], positions[(minID + 1):(maxID)], movingAverageProportion2[(minID + 
                1):(maxID)], col = colramp2[colorId], lty = 1, 
                lwd = 2)
        }
        #legend("topright", bty = "n", col = c(cond1Color, cond2Color), 
        #    legend = conditionsNames, pch = 16, lty = 0, lwd = 2, horiz = TRUE, cex = 0.8)
    }
    if (numberOfDMRs > 0) {
        for (i in 1:numberOfDMRs) {
            range <- DMRs[[i]][queryHits(findOverlaps(DMRs[[i]], 
                region))]
            if (length(range) > 0) {
                start(range) <- pmax(start(range), start(region))
                end(range) <- pmin(end(range), end(region))
                if (length(which(range$regionType == "gain")) > 
                  0) {
                  rect(start(range)[range$regionType == "gain"], 
                    1.075 + (length(DMRs) - i) * 0.1, end(range)[range$regionType == 
                      "gain"], 1.125 + (length(DMRs) - i) * 0.1, 
                    col = DMRsColor[i], border = DMRsColor[i])
                }
                if (length(which(range$regionType == "loss")) > 
                  0) {
                  rect(start(range)[range$regionType == "loss"], 
                    1.075 + (length(DMRs) - i) * 0.1, end(range)[range$regionType == 
                      "loss"], 1.125 + (length(DMRs) - i) * 0.1, 
                    col = DMRsColor[i], border = DMRsColor[i], 
                    density = 30)
                }
                if (length(range) > 0 & (is.null(range$regionType))) {
                  rect(start(range), (1.075 + (length(DMRs) - 
                    i) * 0.1), end(range), (1.125 + (length(DMRs) - 
                    i) * 0.1), col = DMRsColor[i], border = DMRsColor[i])
                }
            }
            ### DMRs in-gragh context labels
            #text(par("usr")[1], 1.1 + (length(DMRs) - i) * 0.1,
            #  pos = 4, names(DMRs)[i], cex = 0.7
            #)

        }
    }

    if (!is.null(gff)) {
        if (is(gff, "GRanges")) {
            .plotGeneticElements(gff, region, c(geneColor, TEColor))
        }
    }
    if (!is.null(labels)) {
        if (length(labels) >= 1) {
            mtext(labels[1], line = 0.7, adj = 0, cex = 1)
        }
    }

########################################
#### add main title
# text((par("usr")[1] + par("usr")[2]) / 2, 1.45, pos = NULL, main, cex = 0.75, font = 2)
try({
    text(par("usr")[1], 1.45, pos = 4, main, cex = 0.75, font = 2)
})
      
      #### DMRs legend
      #text(par("usr")[1], 1.45, pos = 4, "DMRs", cex = 0.7, font = 2)

      #### add legend
      leg_pos1 = end(region) - (end(region) - start(region)) / 10
      leg_pos2 = end(region) + (end(region) - start(region)) / 75
      text(leg_pos1, 1.45, pos = 2, conditionsNames[1], cex = 0.75, font = 2, col = cond1Color)
      text(leg_pos2, 1.45, pos = 2, conditionsNames[2], cex = 0.75, font = 2, col = cond2Color)
      #text(par("usr")[2], 1.45, pos = 2, conditionsNames[1], cex = 0.75, font = 2, col = cond1Color)
      #text(par("usr")[2], 1.325, pos = 2, conditionsNames[2], cex = 0.75, font = 2, col = cond2Color)

      #### add line
      #lines(c(
      #  par("usr")[1], par("usr")[2]
      #  #par("usr")[1] - (par("usr")[2] - par("usr")[1]) / 100,
      #  #par("usr")[2] + (par("usr")[2] - par("usr")[1]) / 100
      # ),
      # y = c(1,1),
      #  lty = 1,
      #  lwd = 0.75,
      #  col = "gray75"
      #)

########################################

    invisible(NULL)
}

###

#### legend
if (T) {
   #svg("/home/yoyerush/yo/genePlot_101124/legend.svg", width = 2.75, height = 1.5, family = "serif")
x.cord = 0.25
y.cord = 0.9

svg("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/genePlot_101124/legend.svg", width = 2.75, height = 2, family = "serif")

par(mar = c(0, 0, 0, 0))
plot.new()

text(0.175, 0.9, "Hyper", cex = 0.9, font = 2, srt = 25)
text(0.35, 0.9, "Hypo", cex = 0.9, font = 2, srt = 25)

first_4_TE = function(x) {
  legend(x.cord - x, y.cord,
    legend = rep("", 6),
    fill = c("white", "white", "white", "white", "white", "#afad43"),
    border = c("white", "white", "white", "white", "white"),
    density = c(NA, NA, NA, NA, NA, 30), angle = c(NA, NA, NA, NA, NA, 30),
    bty = "n", cex = 1.1
  )
}
first_4_TE(0.045)
first_4_TE(0.09)
first_4_TE(0.135)


legend(x.cord, y.cord,
  legend = c("CG", "CHG", "CHH", "", "CDS", "TEs"),
  fill = c("#009E73", "#0072B2", "#CC79A7", "white", "#000000", "#afad43"),
  border = c("#009E73", "#0072B2", "#CC79A7", "white","#000000", "white"),
  density = c(30, 30, 30, NA, NA, 30), angle = c(30, 30, 30, NA, NA, 30),
  bty = "n", cex = 1.1, x.intersp = 1.5
)

legend(x.cord-0.175, y.cord,
  legend = rep("", 6),
  fill = c("#009E73", "#0072B2", "#CC79A7", "white", "#000000", "#afad43"),
  border = c("#009E73", "#0072B2", "#CC79A7", "white", "#000000", "white"),
  density = c(NA, NA, NA, NA, NA, 30), angle = c(NA, NA, NA, NA, NA, 30),
  bty = "n", cex = 1.1
)

# line for CDS
lines(c(0.15, 0.34), c(0.3, 0.3), lty = 1, lwd = 1.5, col = "#000000")

# add box to 'TEs'
rect(0.14, 0.1625, 0.3675, 0.215, border = "#afad43", lwd = 1.5)

dev.off()

}




########################## old
if (F) {
  # svg("/home/yoyerush/yo/genePlot_101124/legend.svg", width = 2.75, height = 1.5, family = "serif")
  svg("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/genePlot_101124/legend.svg", width = 2.75, height = 1.5, family = "serif")
  par(mar = c(0, 0, 0, 0))
  plot.new()

  legend("center",
    legend = c("Hyper-DMRs", "Hypo-DMRs", "Genes", "TEs"),
    fill = c("#009E73", "#009E73", "#000000", "#afad43"),
    border = c("#009E73", "#009E73", "white", "white"),
    density = c(NA, 30, NA, 30), angle = c(NA, 30, NA, 30),
    bty = "n", cex = 1.1, x.intersp = 0.5
  )

  legend("center",
    legend = rep("", 4),
    fill = c("#0072B2", "#0072B2", "#000000", "#afad43"),
    border = c("#0072B2", "#0072B2", "white", "white"),
    density = c(NA, 30, NA, 30), angle = c(NA, 30, NA, 30),
    bty = "n", cex = 1.1, x.intersp = 8.5
  )

  legend("center",
    legend = rep("", 4),
    fill = c("#CC79A7", "#CC79A7", "#000000", "#afad43"),
    border = c("#CC79A7", "#CC79A7", "white", "white"),
    density = c(NA, 30, NA, 30), angle = c(NA, 30, NA, 30),
    bty = "n", cex = 1.1, x.intersp = 10.35
  )

  # add box to the 'Genes' and 'TEs' legend
  rect(0.155, 0.3825, 0.3275, 0.452, border = "#000000", lwd = 1.5)
  rect(0.155, 0.225, 0.3275, 0.3, border = "#afad43", lwd = 1.5)
  lines(c(0.155, 0.3275), c(0.41725, 0.41725), lty = 1, lwd = 1.5, col = "#000000")

  dev.off()
}