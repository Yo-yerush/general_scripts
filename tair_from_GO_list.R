get_tair_from_go <- function(
  go_ids,
  include_offspring = FALSE, ontology = c("BP", "MF", "CC")
) {
    ontology <- match.arg(ontology)

    if (include_offspring) {
        offspring_db <- switch(ontology,
            BP = GO.db::GOBPOFFSPRING,
            MF = GO.db::GOMFOFFSPRING,
            CC = GO.db::GOCCOFFSPRING
        )
        go_ids <- unique(c(go_ids, unlist(lapply(go_ids, function(go_id) as.character(offspring_db[[go_id]])))))
        go_ids <- go_ids[!is.na(go_ids)]
    }
    res <- AnnotationDbi::select(org.At.tair.db::org.At.tair.db, keys = unique(go_ids), keytype = "GO", columns = c("TAIR", "GO", "SYMBOL", "GENENAME"))

    unique(na.omit(res$TAIR))
}