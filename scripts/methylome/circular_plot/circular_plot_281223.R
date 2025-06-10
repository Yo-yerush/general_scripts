
library(circlize)

# Your example dataframe
df <- Chr_df
df_pos = (df$start+df$end)/2
# Read your GFF3 file (modify this to the correct path and format)
#gff3_data <- as.data.frame(genes)[,1:3]
#gff3_pos = (gff3_data$start+gff3_data$end)/2

circos.genomicInitialize(as.data.frame(gff3.trimmed)[,1:3])

# First track - meth_log2FC
circos.trackPlotRegion(df$seqnames, df$start, df$end, ylim = c(min(df$meth_log2FC),max(df$meth_log2FC)), panel.fun = function(region, value, ...) {
  circos.points(df_pos, df$meth_log2FC, col = "#FF0000", pch = ".") # Adjust color and line width as needed
  circos.lines(df_pos, rep(0, length(df$meth_log2FC)), col = "black", pch = 20)
})

# Second track - RNA_log2FC
circos.genomicDensity(as.data.frame(genes)[1:3], count_by = "number")

circos.trackPlotRegion(df$seqnames, df$start, df$end, ylim = c(min(df$RNA_log2FC),max(df$RNA_log2FC)), panel.fun = function(region, value, ...) {
  circos.points(df_pos, df$RNA_log2FC, col = "#0000FF", pch = ".") # Adjust color and line width as needed
  circos.lines(df_pos, rep(0, length(df$RNA_log2FC)), col = "black", pch = 1)
})

# Third track - GFF3 data
circos.genomicDensity(as.data.frame(genes)[1:3], count_by = "number")

# Finish the plot
circos.clear()



circos.genomicRainfall(list(meth_file[meth_file$meth_log2FC > 0,2:4],
                            meth_file[meth_file$meth_log2FC < 0,2:4]),
                       pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))
circos.genomicDensity(meth_file[meth_file$meth_log2FC > 0,2:4], count_by = "number")
circos.genomicDensity(meth_file[meth_file$meth_log2FC < 0,2:4], count_by = "number")
