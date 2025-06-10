library(XML)
library(dplyr)

# Function to extract 'name', 'x', and 'y' from KGML file and return a data frame
extract_data_from_KGML <- function(file_path) {
  # Load XML file
  xml_data <- xmlParse(file_path)

  # Extract <entry> nodes
  entries <- getNodeSet(xml_data, "//entry")
  
  names <- c()
  x_vals <- c()
  y_vals <- c()
  
  i=1
  while (i <= length(entries)) {
    
    # Extract 'name' attribute
    entry_name <- xmlGetAttr(entries[[i]], "name")
    
    # Extract corresponding <graphics> node for 'x' and 'y'
    graphics_node <- getNodeSet(entries[[i]], "./graphics")[1] # Assuming there's only one <graphics> node per <entry>
    x <- xmlGetAttr(graphics_node[[1]], "x")
    y <- xmlGetAttr(graphics_node[[1]], "y")
    
    names <- c(names, entry_name)
    x_vals <- c(x_vals, x)
    y_vals <- c(y_vals, y)
    
    i=i+1
  }
  
  data_frame <- data.frame(name = names, x = as.numeric(x_vals), y = as.numeric(y_vals))
  return(data_frame)
}

file_path <- "P:/yonatan/ath00270 (1).xml"
data <- extract_data_from_KGML(file_path)
print(data)
