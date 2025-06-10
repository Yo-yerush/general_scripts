library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(stringr)
library(multcompView)
library(rlang)
EV = c("EV 2.1","EV 2.2","EV 3.1","EV 3.2","EV 3.3")
bHLH94 = c("bHLH94 1.1","bHLH94 1.2","bHLH94 2","bHLH94 3.2","bHLH94 4")




path = setwd("C:/Users/yonye/OneDrive/שולחן העבודה/מעבדה/RNAseq_yonatan_2021/DESeq2/transgenic/statistics/bHLH94_VS_EV")

table = read.csv("norm.bHLH94_VS_EV.DE.csv")

gene.name = "XM_031526257.1"

title_main = "gallate 1-beta-glucosyltransferase-like"

comp = c(EV, bHLH94)



title_x = ""  # or # names(data[1])
title_y = "Relative Expression (normalization unit)"  # or # names(data[2])




data = data.frame(colnames(table), unlist(table[table$X == gene.name,], use.names = F))[-1,]
names(data) = data.frame(colnames(table), unlist(table[table$X == gene.name,], use.names = F))[1,]
data[,2] = as.numeric(data[,2])
data$X = comp


X = names(data)[1]
Y = names(data[2])
level.order = unique(data$X)


#########################
data %>%
  ggplot(aes(x=factor(x = !!sym(X), level=level.order), y = !!sym(Y))) +
  geom_bar(stat = "summary",color="black", fun.y = "mean", fill="grey70", width=0.5) +
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
  scale_y_continuous(expand = c(0,0), limits = c(0,20000))


