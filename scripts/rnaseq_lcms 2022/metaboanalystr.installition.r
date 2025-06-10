#install.packages("pacman")

#library(pacman)

#pacman::p_load(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea"))

devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)

#install.packages("ellipse")
#install.packages("rjson")

library("MetaboAnalystR")