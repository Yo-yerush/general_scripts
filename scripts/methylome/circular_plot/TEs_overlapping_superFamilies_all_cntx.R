library(dplyr)
library(RColorBrewer)
library(rtracklayer)
library(circlize)

treatment = "mto1"

### GFF RNAseq, and TEs files
{
  trimm_Chr <- function(gr_obj) {
    remove_seqnames = c("ChrM","ChrC")
    gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
    seqlevels(gr_obj) <- setdiff(seqlevels(gr_obj), remove_seqnames)
    return(sort(gr_obj))
  }
  
  ############# gff3 file #############
  gff3 = import.gff3("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 gff3/TAIR10_GFF3_genes.gff")
  gff3.trimmed = trimm_Chr(gff3)
  max_pos_chr5 = gff3.trimmed[which(gff3.trimmed$type == "chromosome")]@ranges@width[5]
  #genes = gff3.trimmed[which(gff3.trimmed$type == "gene")]
  
  #############  RNAseq res file   #############
  RNA_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/",treatment,"_vs_wt/all.transcripts.",treatment,"_vs_wt.DE.csv")) %>%
    filter(padj < 0.05) %>%
    select(locus_tag)
  transcripts = gff3.trimmed[which(gff3.trimmed$type == "mRNA")] %>%
    as.data.frame() %>%
    mutate(locus_tag = as.character(Parent)) %>%
    select(seqnames,start,end,locus_tag) %>%
    distinct(locus_tag, .keep_all = T)
  RNA_transcriptsLoc = merge.data.frame(transcripts, RNA_file, by = "locus_tag")[,c(2:4,1)]
  
  ############# TE file #############
  TE = read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 transposable elements/TAIR10_Transposable_Elements.txt", sep = "\t")
  TE_4_dens = TE[,c("Transposon_Name","Transposon_min_Start","Transposon_max_End","Transposon_Super_Family")]
  TE_4_dens$Transposon_ID = TE_4_dens$Transposon_Name
  for (chr.i in 1:5) {
    TE_4_dens$Transposon_Name = gsub(paste0("AT",chr.i,"TE.*"), paste0("Chr",chr.i),TE_4_dens$Transposon_Name)
  }
  
  ############# TAIR9 gaps #############
  tair9_gaps = rtracklayer::import.gff("C:/Users/yonye/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR9/tair9_Assembly_gaps.gff")
  tair9_gaps = resize(tair9_gaps, width = width(tair9_gaps) + 100000, fix = "center") %>%
    reduce() %>%  as.data.frame() %>% .[,1:3]
  
  ############# Heterochromatin positions #############
  heteroChr = data.frame(Chr = paste0("Chr",c(1:5)),
                         start = c(12500000, 1250000, 11000000, 1666667, 9444444),
                         end = c(17500000, 7500000, 16250000, 7000000, 15000000))
  
  ############# Centromere positions #############
  cenChr = data.frame(Chr = paste0("Chr",c(1:5)),
                      start = c(14476796, 3462971, 13780083, 3177188, 11207348),
                      end = c(15081019, 3650512, 14388500, 3248799, 11555278))
  #cenChr = data.frame(Chr = paste0("Chr",c(1:3,3,3:5)),
  #                    start = c(15086046, 3607930, 13799418, 13587787, 14208953, 3956022, 11725025),
  #                    end = c(15087045, 3608929, 13800417, 13588786, 14209952, 3957021, 11726024))
  #cenChr = data.frame(Chr = paste0("Chr",c(1:5)),
  #                    start = c(15086046, 3607930, 13799418, 3956022, 11725025),
  #                    end = c(15087045, 3608929, 13800417, 3957021, 11726024))
  #cenChr$start = cenChr$start - 4.9e4
  #cenChr$end = cenChr$end + 4.9e4
}

cntx_file <- function(context) {
  ############# DMRs file
  dmrs_file = read.csv(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/",treatment,"_vs_wt/DMRs_",context,"_",treatment,"_vs_wt.csv"))
  #dmrs_file = dmrs_file[,c("seqnames","start","end","log2FC")]
  #names(dmrs_file)[4] = "meth_log2FC"
  return(dmrs_file)
}

CG_file = cntx_file("CG")
CHG_file = cntx_file("CHG")
CHH_file = cntx_file("CHH")

# bind this contexts data frames and make it as GenomicRanges object. each position make it 1 bp
CG_gr = makeGRangesFromDataFrame(CG_file, keep.extra.columns = T)
CHG_gr = makeGRangesFromDataFrame(CHG_file, keep.extra.columns = T)
CHH_gr = makeGRangesFromDataFrame(CHH_file, keep.extra.columns = T)

#CG_gr = resize(CG_gr, width = 1)
#CHG_gr = resize(CHG_gr, width = 1)
#CHH_gr = resize(CHH_gr, width = 1)

all_cntx_gr = c(CG_gr, CHG_gr, CHH_gr)
all_cntx_gr = all_cntx_gr[!duplicated(all_cntx_gr)]

####### TE and DMRs overlay #######
TE_gr = TE_4_dens
names(TE_gr)[1:3] = c("seqnames","start","end", "Transposon_ID")
TE_gr$strand = "*"
#TE_gr = TE_gr[which(TE_gr$Transposon_Super_Family %in% c("LTR/Copia","LTR/Gypsy","LINE/L1","DNA/MuDR","RC/Helitron")),]
TE_gr$col = "#00000090"
#TE_gr$col[TE_gr$Transposon_Super_Family == "LTR/Copia"] = "#17870e90"
#TE_gr$col[TE_gr$Transposon_Super_Family == "LINE/L1"] = "#ab640290"
TE_gr = makeGRangesFromDataFrame(TE_gr, keep.extra.columns = T)
m = findOverlaps(TE_gr, all_cntx_gr)
overlap_TEs = TE_gr[queryHits(m)]
overlap_TEs = resize(overlap_TEs, width = 1)
overlap_TEs = overlap_TEs[!duplicated(overlap_TEs)]
overlap_TEs_df = as.data.frame(overlap_TEs) %>% distinct(Transposon_ID, .keep_all = T)
########## test
#overlap_TEs_df = rbind(read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/Transposable_Elements_CG_genom_annotations.csv"),
#                       read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/Transposable_Elements_CHG_genom_annotations.csv"),
#                       read.csv("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/Transposable_Elements_CHH_genom_annotations.csv"))
#overlap_TEs_df = overlap_TEs_df %>% distinct(Transposon_Name, .keep_all = T)
##########
family_res_list = list(
  SINE = overlap_TEs_df[grep("^SINE|RathE", overlap_TEs_df$Transposon_Super_Family), 1:3],
  LINE = overlap_TEs_df[grep("^LINE", overlap_TEs_df$Transposon_Super_Family), 1:3],
  Copia = overlap_TEs_df[grep("LTR/Copia", overlap_TEs_df$Transposon_Super_Family), 1:3],
  Helitron = overlap_TEs_df[grep("RC/Helitron", overlap_TEs_df$Transposon_Super_Family), 1:3],
  TIR = overlap_TEs_df[grep("^DNA", overlap_TEs_df$Transposon_Super_Family), 1:3],
  Gypsy = overlap_TEs_df[grep("LTR/Gypsy", overlap_TEs_df$Transposon_Super_Family), 1:3]
)
### Calculate densities for RNAseq and TEs
density_RNA <- genomicDensity(RNA_transcriptsLoc[,1:3], window.size = 1e6, count_by = "number")
density_TE <- genomicDensity(TE_4_dens[,1:3], window.size = 1e6, count_by = "number")
# new scaling
density_RNA$value = density_RNA$value / max(density_RNA$value) * 0.4
density_TE$value = density_TE$value / max(density_TE$value)
dir.create("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/TEs_overlapping_superFamilies_all_cntx/", showWarnings = F)
dir.create(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/TEs_overlapping_superFamilies_all_cntx/",treatment), showWarnings = F)
#####################################
############# the plot #############
svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/TEs_overlapping_superFamilies_all_cntx/",
           treatment,"/TEs_overlapping_superFamilies_DMRs_",treatment,".svg"),
    width = 4.25, height = 4.25, family = "serif")
circos.par(gap.degree = c(rep(4,4),35), start.degree = 90)
circos.genomicInitialize(as.data.frame(gff3.trimmed)[,1:3], sector.names = paste0("Chr ",1:5), axis.labels.cex = 0.4, labels.cex = 1.35)
for (family.i in 1:length(family_res_list)) {
  density_data <- genomicDensity(family_res_list[[family.i]], window.size = 1e5, count_by = "number")
  ylims <- range(density_data$value) * 1.435
  circos.genomicTrackPlotRegion(density_data, ylim = range(density_data$value), bg.border = NA, track.height = 0.085, track.margin = c(0, 0), panel.fun = function(region, value, ...) {
    chr.n = gsub("Chr","",get.cell.meta.data("sector.index"))
    ### heterocromatin
    circos.rect(heteroChr[chr.n,2], ylims[1], heteroChr[chr.n,3], ylims[2], # xleft, ybottom, xright, ytop
                col = "#fcba0320",
                border = NA
    )
    ### centromere
    circos.rect(cenChr[chr.n,2], ylims[1], cenChr[chr.n,3], ylims[2], # xleft, ybottom, xright, ytop
                col = "#fcba0360",
                border = NA
    )
    ### gaps
    #gaps.chr <- tair9_gaps[tair9_gaps$seqnames == paste0("Chr", chr.n), ]
    #for (gap.i in nrow(gaps.chr)) {
    #  circos.rect(gaps.chr[gap.i,2], ylims[1], gaps.chr[gap.i,3], ylims[2], # xleft, ybottom, xright, ytop
    #              col = "#0ca81e90",
    #              border = NA
    #  )
    #}
    #### knob1
    #if (chr.n == 1) {
    #  circos.rect(24000000, ylims[1], 24500000, ylims[2], # xleft, ybottom, xright, ytop
    #              col = "#f5425120",
    #              border = NA
    #  )
    #}
    ### density lines
    colors <- ifelse(value > 5, "red3", "gray15")
    circos.genomicLines(region, value, col = colors,
                        border = TRUE, lty = 1, lwd = 0.5, type = "h")
  })
  ### y-axis labels
  circos.text("Chr1", x = 0, y = 0.5, labels = paste0(names(family_res_list)[family.i],"  "), facing = "downward", cex = 0.6, adj = c(1,0))
}
### DEGs and TEs
# DEGs
circos.genomicTrackPlotRegion(density_RNA, ylim = c(0,1), bg.border = NA, track.height = 0.2, track.margin = c(0, 0), panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, col = "gray80", border = TRUE, type = "l", area = T)})
# TEs
circos.genomicTrackPlotRegion(density_TE, bg.border = NA, track.height = 0.05, track.margin = c(0, 0), panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, col = "#fcba0320", border = TRUE, type = "l", area = T,
                      track.index = get.cell.meta.data("track.index")-1 )})
#circos.genomicDensity(list(RNA_transcriptsLoc[,1:3], TE_4_dens[,1:3]), bg.col = "white", # "#fafcff",
#                      bg.border = NA, count_by = c("percent", "number"), ylim.force = T,
#                      col =c("gray80","#fcba0320"), border = T, track.height = 0.2, track.margin=c(0, 0))
circos.clear()
dev.off()
############################################################################################################
#half_val = max(density_data$value) / 2
#q3_value = half_val + half_val/2  #quantile(density_data$value, 0.75)
#colors <- ifelse(density_data$value > 5, "red", "gray40")

#cols.freq =  genomicDensity(family_res_list[[family.i]], count_by = "number", window.size = 1e5)
#circos.genomicDensity(family_res_list[[family.i]], bg.col = "white", # "#fafcff", 
#                      bg.border = NA, count_by = "number", window.size = 1e5,
#                      col =rep("black", nrow(cols.freq)), border = T, track.height = 0.085, track.margin=c(0, 0), type ="h") #("l", "o", "h", "s")

#circos.text("Chr1", x = 0, y = 0.5, labels = "TEs", facing = "downward", adj = c(1.2, -2.25), cex = 0.6)


### vertical line for heterochromatin rages
#for (het in 1:nrow(heteroChr)) {
#  for (t.indx in 2:(get.cell.meta.data("track.index")-1)) {
#    ylims <- get.cell.meta.data("ylim", sector.index = heteroChr[het,1], track.index = t.indx) * 1.5
#    # Drawing the background rectangle
#    circos.rect(heteroChr[het,2],
#                ylims[1],
#                heteroChr[het,3],
#                ylims[2], # xleft, ybottom, xright, ytop,
#                sector.index = heteroChr[het,1],
#                track.index = t.indx,
#                col = "#4287f510",
#                border = NA
#                )
#    # Drawing the lines
#    for (strt.end in 2:3) {
#      circos.lines(
#        x = rep(heteroChr[het, strt.end], 2),
#        y = ylims,
#        sector.index = heteroChr[het,1],
#        track.index = t.indx,
#        col = "blue4",
#        lwd = 0.5,
#        lty = 2
#      )
#    }
#    
#  }
#}
############################################################################################################

# legends
if (F) {
  counts_norm_hyper = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  counts_norm_hypo = hist(rnorm(1000, mean=0, sd=1), plot = F, breaks = 30)$counts
  rndm_norm_hyper = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0,10))
  rndm_norm_hypo = c(rep(0,10), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
  
  #rndm_norm_TE = c(rep(0,10),(counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), rep(0,10))
  rndm_norm_TE = c((counts_norm_hyper - min(counts_norm_hyper)) / (max(counts_norm_hyper) - min(counts_norm_hyper)), (counts_norm_hypo - min(counts_norm_hypo)) / (max(counts_norm_hypo) - min(counts_norm_hypo)))
  #rndm_norm_TE = c(rep(0,10), rndm_norm_hypo[-c(1:10, (length(rndm_norm_hypo)-10):length(rndm_norm_hypo))], rep(0,10))
  
  rndm_norm_genes = (1-rndm_norm_TE)/1.4
  #rndm_norm_genes = 0.6-rndm_norm_TE
  #rndm_norm_genes[rndm_norm_genes<0.12] = 0.12
  
  rndm = runif(12, min = 0, max = 1)
  legend_titles_circular = c("DMRs","  Group sig. genes","Methylation values","Expression values","  Overlapping TEs","Genome annotation")
  tracks_total_size = 16
  legend_tracks_pos = as.character((tracks_total_size/4)+1)
  text_tracks_pos = as.character((tracks_total_size/4))
  
  colfunc <- colorRampPalette(c("red", "white", "blue"))
  colfunc_vec_1 = colfunc(16)[-c((round(16/2)-1):(16/2), ((16/2)+1):((16/2)+1+1))]
  colfunc_vec_2 = colfunc(100)[-c((round(100/2)-5):(100/2)-1, ((100/2)+1):((100/2)+1+5))]
  
  svg(paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/circular_plot_res/DMRs_all_cntx_TEs/legends.svg"), width = 6, height = 5.75, family = "serif")
  #par(mfrow=c(3,1))
  layout(matrix(1:5,ncol=1,nrow=5),height = c(1,0.5,0.15,0.15,0.15))
  
  par(mar = c(0.1, 0.1, 0.1, 0.1) )
  circos.par("track.height" = 0.7, "canvas.xlim" = c(0, 1), "canvas.ylim" = c(0, 1), "gap.degree" = 0, "clock.wise" = FALSE)
  circos.initialize(factors = as.character(1:tracks_total_size), xlim = c(0, 1)) 
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, area = T, border = T, col = "#FF000080")
  circos.lines(seq(0,1,1/(length(rndm_norm_hypo)-1)), rndm_norm_hypo, area = T, border = T, col = "#304ed180")
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  circos.text(text_tracks_pos, x = 0.1, y = 0.5, labels = "CG context", facing = "downward", ad = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, area = T, border = T, col = "#FF000080")
  circos.lines(seq(0,1,1/(length(rndm_norm_hypo)-1)), rndm_norm_hypo, area = T, border = T, col = "#304ed180")
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  circos.text(text_tracks_pos, x = 0.1, y = 0.5, labels = "CHG context", facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.175) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, area = T, border = T, col = "#FF000080")
  circos.lines(seq(0,1,1/(length(rndm_norm_hypo)-1)), rndm_norm_hypo, area = T, border = T, col = "#304ed180")
  circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  circos.text(text_tracks_pos, x = 0.1, y = 0.5, labels = "CHH context", facing = "downward", adj = c(0, 0.5))
  
  circos.trackPlotRegion(factors = 1:tracks_total_size, ylim = c(0, 1), bg.border = NA, track.height = 0.15) 
  circos.updatePlotRegion(sector.index = legend_tracks_pos, bg.border = NA , bg.col = "#fafcfc", bg.lwd=0.2)
  circos.lines(seq(0,1,1/(length(rndm_norm_TE)-1)), rndm_norm_TE, area = T, border = T, col = "#fcba0330")
  #circos.lines(seq(0,1,1/(length(rndm_norm_hyper)-1)), rndm_norm_hyper, border = T)
  circos.text(text_tracks_pos, x = 0.1, y = 0.5, labels = "TEs density", facing = "downward", adj = c(0, 0.5))
  
  # DMRs
  legend(-0.09,0.125,
         legend=c(substitute(paste(bold("Hyper"),"-DMRs")),
                  substitute(paste(bold("Hypo"),"-DMRs")),
                  "Shared"),
         fill = c("#FF000095","#304ed195","#5e0d3d99"),
         bty="n")
  
  circos.clear()
  dev.off()
}
