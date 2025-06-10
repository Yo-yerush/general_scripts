peel.gene.exp = function(transcript.id, path.to.all_high_VS_low.folder = "C:/Users/yonye/OneDrive/שולחן העבודה/מעבדה/RNAseq_yonatan_2021/DESeq2/peel/statistics/all_high_VS_low/")
{
library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(stringr)
library(multcompView)
library(rlang)


path = setwd(path.to.all_high_VS_low.folder)

table = read.csv("norm.all_high_VS_low.DE.csv")
col.data = read.table("../../coldata.peel.txt", header = T)
col.data$sample = as.character(col.data$sample)
table.names = read.csv("../../../transcript_id.and.description.with.GO.csv")


###
###

gene.name = transcript.id

###
###

title_x = ""  # or # names(data[1])
title_y = "Relative Expression (normalization unit)"  # or # names(data[2])


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
data$X = col.data$sample


X = names(data)[1]
Y = names(data[2])
level.order = unique(data$X)

###############

stat = aov( as.formula(paste( Y,'~',X)), data) 
tukey = TukeyHSD(stat, X, conf.level=0.95) 

Tukey.letters = function(tukey, x.ax){
  Tukey.levels <- tukey[[x.ax]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels$X=rownames(Tukey.labels)
  Tukey.labels
}
LABELS.fun = Tukey.letters(tukey, X) 
LABELS=LABELS.fun[match(level.order, LABELS.fun$X),] 

#########################

labels.mean = data %>% group_by(!!sym(X)) %>% summarise_at(vars(!!sym(Y)), list(name = mean))
labels.sd = data %>% group_by(!!sym(X)) %>% summarise_at(vars(!!sym(Y)), list(name = sd))
level.mean.order = labels.mean[match(level.order, labels.mean[[1]]),] 
level.sd.order = labels.sd[match(level.order, labels.sd[[1]]),] 

y.labels.order = level.mean.order$name+level.sd.order$name

#########################
plot.fun = data %>%
  ggplot(aes(x=factor(x = !!sym(X), level=level.order), y = !!sym(Y))) +
  geom_bar(stat = "summary",color="black" , fun.y = "mean", fill="grey70", width=bar.width) +
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
  
  # geom_jitter(color="steelblue4", size=1, alpha=0.9, width = 0.18) +   
  
  
  scale_y_continuous(expand = c(0,0), limits = c(0,max(y.labels.order)*y.axis.hight)) +  ### change y.axis hight -> limits=c(0,hight))
  geom_text(data=LABELS, aes(x= X, y=(y.labels.order*y.label.hight), label=Letters), vjust = 0, fontface = "bold") ### change tukey labels hight - y=(hight)
#plot(title_main)
plot.fun
}