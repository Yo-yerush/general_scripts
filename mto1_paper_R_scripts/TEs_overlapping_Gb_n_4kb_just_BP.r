library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(GenomicRanges)
library(topGO)

dir.create("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_TEs_overlap_Gb_n_4kb/", showWarnings = F)

for (TE_SF_loop in c("all", "Gypsy", "Copia", "LINE", "Helitron", "DNA")) {

    ########################

    ### part 1 ###
    
    ###### TEs results
    TE_1 <- rbind(
        read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CG/Transposable_Elements_CG_genom_annotations.csv"),
        read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHG/Transposable_Elements_CHG_genom_annotations.csv"),
        read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/CHH/Transposable_Elements_CHH_genom_annotations.csv")
    )
        # dplyr::rename(gene_id = Transposon_Name) %>%
        if (TE_SF_loop != "all") {
            TE_1 = TE_1 %>% filter(grepl(TE_SF_loop, Transposon_Super_Family))
        }
        TE_1 = makeGRangesFromDataFrame(TE_1, keep.extra.columns = T) %>%
        sort()


    ################################
    ###### DEGs
    DEGs <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv")
    names(DEGs) <- gsub("locus_tag", "gene_id", names(DEGs))
    symbol_indx <- DEGs %>%
        dplyr::select(gene_id, gene, log2FoldChange, pValue, short_description) %>%
        dplyr::rename(Symbol = gene)
    symbol_indx$log2FoldChange <- round(symbol_indx$log2FoldChange, 3)

    upregulated <- DEGs %>%
        filter(log2FoldChange > 0 & padj < 0.05) %>% # & gene_model_type == "protein_coding") %>%
        dplyr::select(gene_id)

    downregulated <- DEGs %>%
        filter(log2FoldChange < 0 & padj < 0.05) %>% # gene_model_type == "protein_coding") %>%
        dplyr::select(gene_id)

    non_sig <- DEGs %>%
        filter(padj > 0.5 & gene_model_type == "protein_coding") %>%
        dplyr::select(gene_id)


    ################################
    ###### TAIR10 annotations
    gff3 <- rtracklayer::import.gff3("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10/TAIR10 gff3/TAIR10_GFF3_genes.gff")
    gff3_df <- gff3 %>%
        as.data.frame() %>%
        dplyr::rename(gene_id = ID) %>%
        filter(type == "gene") %>%
        dplyr::select(seqnames, start, end, width, strand, gene_id)


    ################################
    ################################
    ###### overlapp with TEs up/down-stream

    # overlap with TEs function
    overlap_fun <- function(x.gr, x.TE = TE_1) {
        m <- findOverlaps(x.gr, x.TE)
        x.out <- x.gr[queryHits(m)] %>%
            as.data.frame() %>%
            mutate(tmp = paste(.$seqnames, .$start, .$end, .$strand, .$gene_id, sep = "_")) %>%
            distinct(tmp, .keep_all = T) %>%
            dplyr::select(gene_id)

        # each gene and its TEs overlapped family
        TE2gene_df <- data.frame(
            gene_id = x.gr[queryHits(m)]$gene_id,
            TE_family = x.TE[subjectHits(m)]$Transposon_Family,
            TE_super_family = x.TE[subjectHits(m)]$Transposon_Super_Family
        ) %>%
            mutate(tmp = paste(.$gene_id, .$TE_family, .$TE_super_family, sep = "_del_")) %>%
            distinct(tmp, .keep_all = T) %>%
            dplyr::select(-tmp)

        return(list(
            overlapIDs = x.out,
            TE2gene = TE2gene_df
        ))
    }

    # remove suplicates fanilies or super-families in columns
    remove_dup <- function(y) {
        y <- as.character(unique(unlist(strsplit(y, "; "))))
        paste(y, collapse = "; ")
    }

    #######
    ## main function
    DEG_near_TEs <- function(deg_type, TE.upstream = F, TE.downstream = F, TE.both_sides = F) {
        deg.gr <- merge(deg_type, gff3_df) %>%
            relocate(gene_id, .after = last_col()) %>%
            makeGRangesFromDataFrame(., keep.extra.columns = T)
        upstream <- shift(deg.gr, 2000)
        downstream <- shift(deg.gr, -2000)

        up_overlap <- overlap_fun(upstream)$overlapIDs
        down_overlap <- overlap_fun(downstream)$overlapIDs

        upstream_TE_fam <- overlap_fun(upstream)$TE2gene
        downstream_TE_fam <- overlap_fun(downstream)$TE2gene
        names(upstream_TE_fam) <- gsub("TE", "upstream", names(upstream_TE_fam))
        names(downstream_TE_fam) <- gsub("TE", "downstream", names(downstream_TE_fam))

        if (TE.both_sides) {
            return_df <- merge(up_overlap, down_overlap) %>%
                merge(., upstream_TE_fam) %>%
                merge(., downstream_TE_fam)
        } else if (TE.upstream) {
            return_df <- up_overlap %>%
                merge(., upstream_TE_fam)
        } else if (TE.downstream) {
            return_df <- down_overlap %>%
                merge(., downstream_TE_fam)
        }

        # group families and super-families by 'gene_id'
        return_df <- return_df %>%
            group_by(gene_id) %>%
            summarise(
                across(
                    contains("stream_family") | contains("stream_super_family"),
                    ~ remove_dup(paste(., collapse = "; "))
                )
            ) %>%
            merge(symbol_indx, ., by = "gene_id") %>%
            as.data.frame() %>%
            arrange(pValue)

        return(return_df)
    }

    xl_list <- list(
        Up.Stream = rbind(
            DEG_near_TEs(upregulated, TE.upstream = T),
            DEG_near_TEs(downregulated, TE.upstream = T)
        ),
        Down.Stream = rbind(
            DEG_near_TEs(upregulated, TE.downstream = T),
            DEG_near_TEs(downregulated, TE.downstream = T)
        )#,
        #Both.Sides = rbind(
        #    DEG_near_TEs(upregulated, TE.both_sides = T),
        #    DEG_near_TEs(downregulated, TE.both_sides = T)
        #)
    )

    ########################

    ########################

    ### part 2 ###

    ########################
    #upstream <- readxl::read_excel("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_near_DEGs/TEs_near_DEGs_201124.xlsx", sheet = "Up.Stream")
    upstream = xl_list[["Up.Stream"]]
    #if (TE_SF_loop != "all") {
    #    upstream <- upstream %>% filter(grepl(TE_SF_loop, upstream_super_family))
    #}
    upstream <- upstream %>% dplyr::select(gene_id)

    #downstream <- readxl::read_excel("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/TEs_near_DEGs/TEs_near_DEGs_201124.xlsx", sheet = "Down.Stream")
    downstream = xl_list[["Down.Stream"]]
    #if (TE_SF_loop != "all") {
    #    downstream <- downstream %>% filter(grepl(TE_SF_loop, downstream_super_family))
    #}
    downstream <- downstream %>% dplyr::select(gene_id)
    ########################

    RNAseq <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/met23/mto1_vs_wt/all.transcripts.mto1_vs_wt.DE.csv") %>%
        dplyr::rename(gene_id = locus_tag, DEG_log2FC = log2FoldChange) %>%
        dplyr::select(gene_id, DEG_log2FC, padj, pValue)

    overlapped_TE_context <- function(context, TE_SF) {
        TE <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Transposable_Elements_", context, "_genom_annotations.csv"))
        if (TE_SF != "all") {
            TE <- TE %>% filter(grepl(TE_SF, Transposon_Super_Family))
        }
        TE <- makeGRangesFromDataFrame(TE, keep.extra.columns = TRUE)

        overlapped <- data.frame()

        gene_feature <- read.csv(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/genome_annotation/", context, "/Genes_", context, "_genom_annotations.csv")) %>%
            makeGRangesFromDataFrame(keep.extra.columns = TRUE)

        m <- findOverlaps(gene_feature, TE)
        m_df <- gene_feature[queryHits(m)]
        m_df$overlapping_TE <- paste(TE$Transposon_Super_Family[subjectHits(m)], TE$Transposon_Family[subjectHits(m)], sep = ", ")
        # m_df$Transposon_Super_Family <- TE$Transposon_Super_Family[subjectHits(m)]
        # m_df$Transposon_Family <- TE$Transposon_Family[subjectHits(m)]

        overlapped <- rbind(overlapped, as.data.frame(m_df)) %>%
            dplyr::select(gene_id)

        return(overlapped)
    }

    ################
    remove_dup_DMR <- function(y) {
        y <- as.character(unique(unlist(strsplit(y, ","))))
        paste(y, collapse = ",")
    }
    ################

    ################
    overlapped_TE_sig_genes <- rbind(
        overlapped_TE_context("CG", TE_SF_loop),
        overlapped_TE_context("CHG", TE_SF_loop),
        overlapped_TE_context("CHH", TE_SF_loop),
        upstream,
        downstream
    ) %>%
        distinct(gene_id) %>%
        merge.data.frame(RNAseq, ., by = "gene_id") %>%
        filter(padj < 0.05) %>%
        mutate(for_geneList = 1) %>%
        dplyr::select(gene_id, for_geneList)


    overlapped_TE <- merge(overlapped_TE_sig_genes, RNAseq, by = "gene_id", all.y = T)

    overlapped_TE$for_geneList[is.na(overlapped_TE$for_geneList)] <- 0

    geneList <- overlapped_TE$for_geneList
    names(geneList) <- overlapped_TE$gene_id


    res_list <- list()
    for (GO_type_loop in c("BP")) { #, "MF", "CC")) {
        myGOdata <- new("topGOdata",
            ontology = GO_type_loop,
            allGenes = geneList,
            geneSelectionFun = function(x) (x == 1),
            # description = "Test",
            annot = annFUN.org,
            # nodeSize = 5,
            mapping = "org.At.tair.db"
        )

        sg <- sigGenes(myGOdata)
        str(sg)
        numSigGenes(myGOdata)

        resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")


        allRes <- GenTable(myGOdata,
            Fisher = resultFisher,
            orderBy = "Fisher", ranksOf = "Fisher", topNodes = length(resultFisher@score)
        )
        allRes$Fisher <- as.numeric(allRes$Fisher)
        allRes$Term <- gsub(",", ";", allRes$Term)
        allRes$type <- GO_type_loop

        res_list[[GO_type_loop]] <- allRes[allRes$Fisher <= 0.01, ]
    }

    ### plot
    # just BP !!!!!!!!!
    res_bind <- rbind(res_list[["BP"]]) # , res_list[["CC"]], res_list[["MF"]])
    res_bind$type <- "Biological Process"

    res_bind <- res_bind[!grepl("cellular_component|biological_process|molecular_function|macromolecule biosynthetic process", res_bind$Term), ]

    # gain_col = "#cf534c"
    # loss_col = "#6397eb"

    bubble_plot <- res_bind %>%
        ggplot(aes(Significant, reorder(Term, Significant), size = Annotated, color = Fisher)) +
        scale_color_gradient("p.value", low = "#F2A672", high = "black") + # theme_classic() +
        labs( # title = paste0(type_name," - ",gain_loss,"regulated transcripts"),
            x = "Significant", y = ""
        ) +
        theme_bw() +
        theme(
            # plot.title=element_text(hjust=0.5),
            # legend.key.size = unit(0.25, "cm"),
            # legend.title = element_text(size = 9.5),
            # legend.position = "right",
            # text = element_text(family = "serif")
            legend.position = "none"
        ) +
        geom_point() +
        facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
        guides(color = guide_colorbar(order = 1, barheight = 4))

    # print(paste(treatment,context,annotation, sep = "_"))

    svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_TEs_overlap_Gb_n_4kb/GO_", TE_SF_loop, "_Gb_n_4kb_DEGs_overlap_with_TEs.svg"), width = 4.25, height = 3, family = "serif")
    print(bubble_plot)
    dev.off()

    svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/GO_TEs_overlap_Gb_n_4kb/GO_", TE_SF_loop, "_Gb_n_4kb_DEGs_overlap_with_TEs_Annotated_legend.svg"), width = 1, height = 1.5, family = "serif")
    legend_plot <- ggplot(res_bind, aes(x = 1, y = 1, size = Annotated)) +
        geom_point() +
        scale_size(
            # range = c(1, 6),
            name = " ",
            breaks = c(min(res_bind$Annotated), 50, 100, max(res_bind$Annotated))
        ) +
        theme_void() +
        theme(legend.position = "right")

    legend_only <- cowplot::get_legend(legend_plot)
    grid.newpage()
    grid.draw(legend_only)
    dev.off()
}
