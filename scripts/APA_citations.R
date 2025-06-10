APA_citations <- function(citation_query) {
  
  library(rcrossref)
  library(RefManageR)
  
  result <- cr_works(query = citation_query)
  
  result$data = result$data[1,]
  
  # Function to format authors
  format_authors <- function(authors) {
    paste(authors$family, substr(authors$given, 1, 1), sep = " ", collapse = "., ")
  }
  #format_authors <- function(authors) {
  #    paste(authors$given, authors$family, sep = " ", collapse = ", ")
  #  }
  
  
  papers <- data.frame(
    title = result$data$title,
    #author = sapply(result$data$author, function(x) x$given[1]),
    author = format_authors(result$data$author[[1]]),
    journal = result$data$container.title,
    year = substr(result$data$published.print, 1, 4),
    volume = result$data$volume,
    issue = result$data$issue,
    pages = result$data$page
    #doi = result$data$DOI
  )
  
  if (is.na(papers$issue) & is.na(papers$pages)) {
    APA_output = paste0(papers$author,". (",papers$year,"). ",papers$title,". ",papers$journal,", ",papers$volume,".")
  } else if (is.na(papers$issue)) {
    APA_output = paste0(papers$author,". (",papers$year,"). ",papers$title,". ",papers$journal,", ",papers$volume,", ",papers$pages,".")
  } else if (is.na(papers$pages)) {
    APA_output = paste0(papers$author,". (",papers$year,"). ",papers$title,". ",papers$journal,", ",papers$volume,"(",papers$issue,").")
  } else {
    APA_output = paste0(papers$author,". (",papers$year,"). ",papers$title,". ",papers$journal,", ",papers$volume,"(",papers$issue,"), ",papers$pages,".")
  }
  
  
  return(APA_output)
  
  
  
  
  
  if (F) {
    citations <- BibEntry(
      bibtype = "Article",
      title = papers$title,
      author = as.person(papers$author),
      #author = papers$author,
      journal = papers$journal,
      year = papers$year,
      volume = papers$volume,
      number = papers$issue,
      pages = papers$pages,
      doi = papers$doi
    )
    apa_citations <- toBiblatex(citations, style = "apa")
    
  }
  
  
  
}