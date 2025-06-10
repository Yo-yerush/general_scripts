create_coloraxis = function(comp_pathways = "",
                            specific_hex = NA,
                            pathway_cor_file = "P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/cor.map/cor.map.3nd/cor_table_metabo/pathway_table/correlation_table.csv",
                            check_color_index = F) {


library(grDevices)

cor.path.table = read.csv(pathway_cor_file)[,1]
cor.path.table = gsub(" X.*","",cor.path.table)
cor.path.table = gsub("lignins", "Monolignols",cor.path.table)
cor.path.table = gsub("leucoanthocyanidin", "Leucoanthocyanidins",cor.path.table)
cor.path.table = gsub("flav", "Flav",cor.path.table)
cor.path.table = gsub("isoFlavones", "Isoflavones",cor.path.table)
cor.path.table = gsub("coumarin", "Coumarins",cor.path.table)
cor.path.table = gsub("HTs", "Hydrolysable Tannins",cor.path.table)
cor.path.table = gsub("proanthocyanidin", "Chatechins",cor.path.table)
cor.path.table = gsub("Anthocyanidin", "Anthocyanins",cor.path.table)

n = length(unique(cor.path.table))
color.hex = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n,
                    alpha = NULL, rev = FALSE)


comp.color.index = data.frame(unique(cor.path.table),color.hex)
#check color.table
if (check_color_index == T) {
  names(comp.color.index)[1] = "pathway.name"
  print(comp.color.index)
}
# specific color for pathway
if (comp_pathways[1] != "" & is.na(specific_hex[1]) == F) {
  names(comp.color.index)[1] = "pathway.name"
  for (i in 1:length(comp_pathways)) {
    comp.color.index[comp.color.index[,1] == comp_pathways[i],][,2] = specific_hex[i]
  }
}
names(comp.color.index) = c("X","Y")

# create table of all metabolites with thire related pathway and color
cor.path.finel.table = data.frame(cor.path.table,cor.path.table)
i=1
while (i <= length(comp.color.index$X)) {
  cor.path.finel.table[,2] = gsub(comp.color.index$X[i],comp.color.index$Y[i],cor.path.finel.table[,2])
  i=i+1
}

############################################
# make only the chose pathway to be colored (all the auther will get white color)
if (comp_pathways[1] != "") {
  u=1
  while (u <= length(cor.path.finel.table$cor.path.table)) {
    if (all(cor.path.finel.table$cor.path.table[u] != comp_pathways) == T) {
      cor.path.finel.table$cor.path.table.1[u] = "#FFFFFF" 
    }
    u=u+1
  }
}

plot.new()
# get the "bar-code" plot
barplot(rep(0.2, length(cor.path.finel.table[,2])), col = cor.path.finel.table[,2], ylim = c(0,2),
        yaxt = "n", border = NA, space = 0, names.arg=1:length(cor.path.finel.table[,2]))
#legend("topleft", legend = comp.color.index$X, pch=16, pt.cex=1.5, cex=0.8, bty='n', col = comp.color.index$Y)
if (comp_pathways[1] == "") {
  legend("top", legend = comp.color.index$X, pch=16, pt.cex=1.5, cex=0.8, bty='n', col = comp.color.index$Y)
} else {
  legend("top", legend = comp_pathways, pch=16, pt.cex=1.5, cex=0.8, bty='n', col = specific_hex)
}
}

###################################

#  #only_legend = F
#if (only_legend == T) {
# get legend plot
#  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#  #legend("topleft", legend = comp.color.index$X, pch=16, pt.cex=1.5, cex=0.8, bty='n', col = comp.color.index$Y)
#  legend("topleft", legend = comp_pathways, pch=16, pt.cex=1.5, cex=0.8, bty='n', col = specific_hex)
#} else {
## get the "bar-code" plot
#  barplot(rep(0.2, length(cor.path.finel.table[,2])), col = cor.path.finel.table[,2],
#          yaxt = "n", border = NA, space = 0, names.arg=1:length(cor.path.finel.table[,2]))
#}

###################################

