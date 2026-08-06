offspring_fun <- function(go_id, xx = as.list(GO.db::GOBPOFFSPRING)) { # 'GOBPCHILDREN' for child terms

  child_terms_0 <- as.character(xx[[go_id]])
  child_terms <- child_terms_0

  for (i in 1:length(child_terms_0)) {
    child_terms <- c(child_terms, as.character(xx[[child_terms[i]]]))
  }

  return(unique(child_terms[!is.na(child_terms)])) # %>% paste(collapse = "|"))
}
