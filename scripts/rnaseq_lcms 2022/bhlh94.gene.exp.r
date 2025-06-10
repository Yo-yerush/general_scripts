bhlh94.gene.exp <- function(transcript.id,
                            path.to.bhlh94_VS_ev.folder = "P:/yonatan/RNAseq_yonatan_2021/DESeq2/transgenic/statistics/bHLH94-like_VS_EV") {
  
library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(stringr)
library(multcompView)
library(rlang)


path = setwd(path.to.bhlh94_VS_ev.folder)

table = read.csv("norm.bHLH94-like_VS_EV.DE.csv")
table_padj = read.csv("transgenic.all.bHLH94-like_VS_EV.DE.csv")[,c(1,4)]
col.data = read.table("../../coldata.transgenic.txt", header = T)
col.data$sample = as.character(col.data$sample)
table.names = read.csv("../../../transcript_id.and.description.with.GO.csv")


###
###

gene.name = transcript.id

###
###


title_x = ""  # or # names(data[1])
title_y = "Relative Expression"  # or # names(data[2])


y.axis.hight = 1.1
y.label.hight = 1

bar.width = 0.5
error.width = 0.1

##############

title_main = paste(unlist(table.names[table.names$transcript_id == gene.name,][,c(1,3)], use.names = F)[2], " (",
                   unlist(table.names[table.names$transcript_id == gene.name,][,c(1,3)], use.names = F)[1], ")",sep = "")

data = data.frame(colnames(table), unlist(table[table$X == gene.name,], use.names = F))[-1,]
names(data) = data.frame(colnames(table), unlist(table[table$X == gene.name,], use.names = F))[1,]
data[,2] = as.numeric(data[,2])
data$X = col.data$sample[c(1:5,16:18,20)]


X = names(data)[1]
Y = names(data[2])
level.order = unique(data$X)

###############

#t.test_p.value = t.test(data[1:5,2],data[6:9,2], var.equal = F)$p.value
t.test_p.value = table_padj[table_padj$transcript_id == gene.name,2]
if (t.test_p.value <= 0.05 & t.test_p.value > 0.01) {
  p_value_lable = "*"
} else if (t.test_p.value <= 0.01 & t.test_p.value > 0.001) {
  p_value_lable = "**"
} else if (t.test_p.value <= 0.001) {
  p_value_lable = "***"
} else if (t.test_p.value > 0.05) {
  p_value_lable = "ns"
}

#########################



#########################
plot.fun = data %>%
  ggplot(aes(x=factor(x = !!sym(X), level=level.order), y = !!sym(Y))) +
  geom_bar(stat = "summary",color="black" , fun.y = "mean", fill=c("#562842", "grey70"), width=bar.width) +
  geom_errorbar (stat = "summary", fun.data = "mean_se", width = error.width) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 10, face="bold"),
        axis.text.y = element_text(size = 10, face="bold"),
        axis.line.x = element_line(size = 1.1),
        axis.line.y = element_line(size = 1.1),
        axis.ticks = element_line(size = 1.1) , 
        axis.ticks.length = unit(0.01, "cm")) +
  labs(title = title_main, x = title_x, y = title_y) +
  
   geom_jitter(color="steelblue4", size=1, alpha=0.9, width = 0.18) +   
  
  
  scale_y_continuous(expand = c(0,0), limits = c(0,max(data[,2]))) +  ### change y.axis hight -> limits=c(0,hight))
  annotate("text", x = 1.5, y = max(data[,2])*0.95, label= p_value_lable, size = 7)
           
  #  stat_compare_means(method = "t.test", label.x = 1.5, label.y = max(data[,2])*0.9, ref.group = "EV")
print(paste("bHLH94-like VS EV: p.value = ", t.test_p.value, sep = ""))
plot.fun
}