# Install required packages if they are not already installed
if (!require("data.tree")) install.packages("data.tree", dependencies = TRUE)
if (!require("DiagrammeR")) install.packages("DiagrammeR", dependencies = TRUE)

# Load libraries
library(data.tree)
library(DiagrammeR)

# Specify the directory path (use double backslashes in Windows paths)
dir_path <- "C:/Users/yonatany/OneDrive - Migal/Desktop/mto1_vs_wt"

# Get all file and directory paths recursively
all_paths <- list.files(path = dir_path, recursive = TRUE, full.names = TRUE)

# Remove paths containing 'genom_annotations.csv' or 'mto1_vs_wt.bedGragh'
all_paths <- all_paths[!grepl("genom_annotations.csv|mto1_vs_wt.bedGragh|topGO.csv|weight01_\\d+_all.pdf|KEGG.csv|GO.svg|features.csv", all_paths)]

# Create relative paths by removing the base directory path
relative_paths <- gsub(dir_path, "mto1_vs_wt", all_paths)

# Build path strings for the data.tree structure
path_strings <- file.path(basename(dir_path), relative_paths, fsep = "/")

# Create the tree structure
tree <- as.Node(data.frame(pathString = path_strings))
tree





# Plot the tree using DiagrammeR
graph <- ToDiagrammeRGraph(tree)
render_graph(graph)
