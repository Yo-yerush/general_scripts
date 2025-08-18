library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(ggbreak)

TE_df <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/Arabidopsis_db/TAIR10_Transposable_Elements.txt", sep = "\t")

TE_df$Transposon_ID <- TE_df$Transposon_Name
TE_df$Transposon_Name <- gsub(paste0("AT1TE.*"), "Chr1", TE_df$Transposon_Name)
TE_df$Transposon_Name <- gsub(paste0("AT2TE.*"), "Chr2", TE_df$Transposon_Name)
TE_df$Transposon_Name <- gsub(paste0("AT3TE.*"), "Chr3", TE_df$Transposon_Name)
TE_df$Transposon_Name <- gsub(paste0("AT4TE.*"), "Chr4", TE_df$Transposon_Name)
TE_df$Transposon_Name <- gsub(paste0("AT5TE.*"), "Chr5", TE_df$Transposon_Name)

# TE_df = TE_df[,-c(5:6)]
names(TE_df)[1:4] <- c("seqnames", "strand", "start", "end")
TE_df$strand <- ifelse(TE_df$strand == "true", "+", "-")
TE_df$strand <- "*"
# TE_df <- TE_df %>% select(-Transposon_ID)# %>% filter(!Transposon_Super_Family == "Unassigned")
TE_gr <- makeGRangesFromDataFrame(TE_df, keep.extra.columns = T) # dont need families characterization here

##########################################

CG <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/DMRs_CG_mto1_vs_wt.csv")[, c("seqnames", "start", "end", "width", "strand", "log2FC")]
CHG <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/DMRs_CHG_mto1_vs_wt.csv")[, c("seqnames", "start", "end", "width", "strand", "log2FC")]
CHH <- read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/DMRs_CHH_mto1_vs_wt.csv")[, c("seqnames", "start", "end", "width", "strand", "log2FC")]

DMRs_gr <- GRangesList(
  CG = makeGRangesFromDataFrame(CG, keep.extra.columns = T),
  CHG = makeGRangesFromDataFrame(CHG, keep.extra.columns = T),
  CHH = makeGRangesFromDataFrame(CHH, keep.extra.columns = T)
)

##########################################

# DMRs overlap with TEs
overlap_fun <- function(te, dmr, is.unique) {
  x <- findOverlaps(te, dmr, ignore.strand = TRUE)
  qh <- queryHits(x)
  sh <- subjectHits(x)
  xx <- te[qh]
  # FIX: assign just the DMR log2FC; also store exact bp overlap for weighting later
  mcols(xx)$log2FC <- mcols(dmr)$log2FC[sh]
  mcols(xx)$size <- width(xx)
  mcols(xx)$iwidth <- width(pintersect(te[qh], dmr[sh]))

  if (is.unique) {
    ## here i only want the average log2FC within TEs. its ok that thee few DMRs on the same TE
    xx <- xx[!duplicated(xx$Transposon_ID)]
    ## if its not unique, we can assume that long TEs will overlap more than short (more DMRs on long TEs cous they are long)
  }
  return(xx)
}

# short VS long TEs
range_widths <- width(TE_gr)

long_TEs <- TE_gr[range_widths > 4000]
short_TEs <- TE_gr[range_widths < 500]

##########################################

overlap_TEs <- GRangesList(
  CG = overlap_fun(TE_gr, DMRs_gr[["CG"]], F),
  CHG = overlap_fun(TE_gr, DMRs_gr[["CHG"]], F),
  CHH = overlap_fun(TE_gr, DMRs_gr[["CHH"]], F)
)

overlap_long <- GRangesList(
  CG = overlap_fun(long_TEs, DMRs_gr[["CG"]], F),
  CHG = overlap_fun(long_TEs, DMRs_gr[["CHG"]], F),
  CHH = overlap_fun(long_TEs, DMRs_gr[["CHH"]], F)
)

overlap_short <- GRangesList(
  CG = overlap_fun(short_TEs, DMRs_gr[["CG"]], F),
  CHG = overlap_fun(short_TEs, DMRs_gr[["CHG"]], F),
  CHH = overlap_fun(short_TEs, DMRs_gr[["CHH"]], F)
)

##########################################

# Calculate the Spearman correlation coefficient
spearman_correlation <- list(
  CG = cor(overlap_TEs[["CG"]]$size, overlap_TEs[["CG"]]$log2FC, method = "spearman"),
  CHG = cor(overlap_TEs[["CHG"]]$size, overlap_TEs[["CHG"]]$log2FC, method = "spearman"),
  CHH = cor(overlap_TEs[["CHH"]]$size, overlap_TEs[["CHH"]]$log2FC, method = "spearman")
)


long_spearman_correlation <- list(
  CG = cor(overlap_long[["CG"]]$size, overlap_long[["CG"]]$log2FC, method = "spearman"),
  CHG = cor(overlap_long[["CHG"]]$size, overlap_long[["CHG"]]$log2FC, method = "spearman"),
  CHH = cor(overlap_long[["CHH"]]$size, overlap_long[["CHH"]]$log2FC, method = "spearman")
)

short_spearman_correlation <- list(
  CG = cor(overlap_short[["CG"]]$size, overlap_short[["CG"]]$log2FC, method = "spearman"),
  CHG = cor(overlap_short[["CHG"]]$size, overlap_short[["CHG"]]$log2FC, method = "spearman"),
  CHH = cor(overlap_short[["CHH"]]$size, overlap_short[["CHH"]]$log2FC, method = "spearman")
)

# average
TE_average <- list(
  CG = mean(overlap_TEs[["CG"]]$log2FC),
  CHG = mean(overlap_TEs[["CHG"]]$log2FC),
  CHH = mean(overlap_TEs[["CHH"]]$log2FC)
)

long_average <- list(
  CG = mean(overlap_long[["CG"]]$log2FC),
  CHG = mean(overlap_long[["CHG"]]$log2FC),
  CHH = mean(overlap_long[["CHH"]]$log2FC)
)

short_average <- list(
  CG = mean(overlap_short[["CG"]]$log2FC),
  CHG = mean(overlap_short[["CHG"]]$log2FC),
  CHH = mean(overlap_short[["CHH"]]$log2FC)
)



##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################



agg_from_overlap <- function(gr, context_label) {
  if (length(gr) == 0L) {
    return(dplyr::tibble())
  }
  df <- as.data.frame(gr)
  # expects columns: Transposon_ID, width (TE length), log2FC, iwidth (added above)
  df %>%
    group_by(Transposon_ID) %>%
    summarise(
      te_len = first(width),
      n_dmr = n(),
      dmr_bp = sum(iwidth, na.rm = TRUE),
      wmean_log2FC = weighted.mean(log2FC, iwidth, na.rm = TRUE),
      hyper_bp = sum(iwidth[log2FC > 0], na.rm = TRUE),
      hypo_bp = sum(iwidth[log2FC < 0], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      context = context_label,
      dmr_cov = dmr_bp / te_len
    )
}

agg_all <- bind_rows(
  agg_from_overlap(overlap_TEs[["CG"]], "CG"),
  agg_from_overlap(overlap_TEs[["CHG"]], "CHG"),
  agg_from_overlap(overlap_TEs[["CHH"]], "CHH")
)

# TE meta and length bins (re-using your thresholds)
te_meta <- tibble(
  Transposon_ID = mcols(TE_gr)$Transposon_ID,
  te_len = width(TE_gr),
  length_bin = case_when(
    te_len > 4000 ~ "long",
    te_len < 500 ~ "short",
    TRUE ~ "mid"
  )
) %>%
  select(-te_len)

agg_all <- agg_all %>%
  left_join(te_meta, by = "Transposon_ID")

##########################################
# NEW: DMR rate per kb by bin and context (useful, length-normalized)
rate_tbl <- agg_all %>%
  filter(length_bin %in% c("long", "short")) %>%
  group_by(length_bin, context) %>%
  summarise(
    dmr_bp_sum = sum(dmr_bp, na.rm = TRUE),
    te_bp_sum = sum(te_len, na.rm = TRUE),
    rate_dmr_per_kb = (dmr_bp_sum / te_bp_sum) * 1000,
    .groups = "drop"
  )

print(rate_tbl)

##########################################
# NEW: Presence/absence enrichment (Fisher) — TEs with ≥1 DMR vs none, long vs short
has_any_tbl <- te_meta %>%
  left_join(agg_all %>% distinct(Transposon_ID, context), by = "Transposon_ID") %>%
  mutate(has_any_dmr = !is.na(context)) %>%
  filter(length_bin %in% c("long", "short"))

fisher_by_context <- function(ctx) {
  tmp <- te_meta %>%
    left_join(agg_all %>% filter(context == ctx) %>% distinct(Transposon_ID) %>% mutate(has_any_dmr = TRUE),
      by = "Transposon_ID"
    ) %>%
    mutate(has_any_dmr = ifelse(is.na(has_any_dmr), FALSE, TRUE)) %>%
    filter(length_bin %in% c("long", "short"))
  tab <- table(tmp$length_bin, tmp$has_any_dmr)
  list(context = ctx, table = tab, fisher = fisher.test(tab))
}
fisher_CG <- fisher_by_context("CG")
fisher_CHG <- fisher_by_context("CHG")
fisher_CHH <- fisher_by_context("CHH")

print(fisher_CG)
print(fisher_CHG)
print(fisher_CHH)

##########################################
# NEW: Length ~ methylation effect (Spearman) using aggregated (per-TE) log2FC
agg_cors <- agg_all %>%
  filter(length_bin %in% c("long", "short")) %>%
  group_by(context) %>%
  summarise(
    rho = suppressWarnings(cor(te_len, wmean_log2FC, method = "spearman", use = "complete.obs")),
    .groups = "drop"
  )

print(agg_cors)

##########################################
# NEW: Quick plot — DMR rate per kb for long vs short by context
ggplot(rate_tbl, aes(x = length_bin, y = rate_dmr_per_kb, fill = context)) +
  geom_col(position = "dodge") +
  labs(x = "TE length bin", y = "DMR bp per kb of TE", title = "DMR density by TE length and context") +
  theme_bw()

##########################################
# (Optional) Edge vs body enrichment for long TEs in CHH (RdDM-at-edges vs CMT2-in-body signal)
# Small, self-contained block; comment out if not needed
edge_body_enrichment <- function(te_set, dmr_gr, pad = 200L) {
  left <- GRanges(seqnames(te_set), IRanges(start(te_set), pmin(end(te_set), start(te_set) + pad - 1)))
  right <- GRanges(seqnames(te_set), IRanges(pmax(start(te_set), end(te_set) - pad + 1), end(te_set)))
  edges <- c(left, right)
  body <- GRanges(seqnames(te_set), IRanges(start(te_set) + pad, pmax(start(te_set) + pad, end(te_set) - pad)))

  he <- findOverlaps(edges, dmr_gr, ignore.strand = TRUE)
  hb <- findOverlaps(body, dmr_gr, ignore.strand = TRUE)

  edge_bp <- sum(width(pintersect(edges[queryHits(he)], dmr_gr[subjectHits(he)])))
  body_bp <- sum(width(pintersect(body[queryHits(hb)], dmr_gr[subjectHits(hb)])))

  edge_len <- sum(width(edges))
  body_len <- sum(width(body))
  tibble(rate_edge = edge_bp / edge_len * 1000, rate_body = body_bp / body_len * 1000)
}

# Example for CHH:
long_idx <- which(width(TE_gr) > 4000)
edge_body_tbl <- edge_body_enrichment(TE_gr[long_idx], DMRs_gr[["CHH"]], pad = 200L)
print(edge_body_tbl)
