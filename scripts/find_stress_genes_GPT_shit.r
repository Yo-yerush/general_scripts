############################################################
# Script: find_stress_genes.R
############################################################

library(org.At.tair.db)
library(GO.db)
library(AnnotationDbi)

############################################################
# 2) Identify all GO terms that contain the word "stress"
############################################################

# Retrieve all GO IDs from GO.db
all_go_ids <- keys(GO.db)

# Extract GO term information: GOID, TERM, and ONTOLOGY
go_info <- select(GO.db,
    keys = all_go_ids,
    columns = c("GOID", "TERM", "ONTOLOGY"),
    keytype = "GOID"
)

# Filter GO terms that include "stress" (case-insensitive)
stress_terms <- go_info[grepl("stress", go_info$TERM, ignore.case = TRUE), ]
head(stress_terms)
# This data frame contains GO IDs and their corresponding names (TERM)
# that mention 'stress'.

############################################################
# 3) Find all Arabidopsis genes annotated to these "stress" GO terms
############################################################

# We will use the GO IDs from stress_terms as keys
stress_go_ids <- unique(stress_terms$GOID)

# Now query org.At.tair.db to find the Arabidopsis genes associated
# with these stress-related GO terms
stress_genes <- select(org.At.tair.db,
    keys = stress_go_ids,
    columns = c("TAIR", "SYMBOL", "GENENAME", "GO"),
    keytype = "GO"
)

# 'stress_genes' will contain rows with:
# - GO: The GO term ID
# - TAIR: TAIR gene ID (e.g. AT1G01010)
# - SYMBOL: Commonly used gene symbol (if available)
# - GENENAME: Full gene name/description
head(stress_genes)

############################################################
# 4) (Optional) Clean up or export the results
############################################################

# Remove duplicates if you only want unique gene entries
stress_genes_unique <- unique(stress_genes[, c("TAIR", "SYMBOL", "GENENAME")])

stress_genes_unique = stress_genes_unique[-1,]
# Examine the first few entries
head(stress_genes_unique)

# You can write this out to a CSV file, for example:
# write.csv(stress_genes_unique, file = "Arabidopsis_stress_genes.csv", row.names = FALSE)

############################################################
# End of Script
############################################################
