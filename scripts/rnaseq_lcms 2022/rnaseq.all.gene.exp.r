rnaseq.all.gene.exp <- function(transcript.id) {

  if (length(dev.list()!=0)) {dev.off()}
source("P:/yonatan/scripts/rnaseq_lcms 2022/peel.gene.exp.r")
source("p:/yonatan/scripts/rnaseq_lcms 2022/clone.gene.exp.r")
source("p:/yonatan/scripts/rnaseq_lcms 2022/bhlh94.gene.exp.r")
  suppressWarnings(
    tryCatch(suppressWarnings(return(list(
      bhlh94.gene.exp(transcript.id),
      clone.gene.exp(transcript.id),
      peel.gene.exp(transcript.id)
    ))), error = function(e) {
        return(list(
          clone.gene.exp(transcript.id),
          peel.gene.exp(transcript.id)
        ))
      }, error = function(e) {
        return(list(
          bhlh94.gene.exp(transcript.id),
          peel.gene.exp(transcript.id)
        ))
      }, error = function(e) {
        return(list(
          clone.gene.exp(transcript.id),
          bhlh94.gene.exp(transcript.id)
        ))
      }, error = function(e) {
        return(list(
          peel.gene.exp(transcript.id)
        ))
      }, error = function(e) {
        return(list(
          clone.gene.exp(transcript.id)
        ))
      }, error = function(e) {
        return(list(
          bhlh94.gene.exp(transcript.id)
        ))
      }
    )
)
}

