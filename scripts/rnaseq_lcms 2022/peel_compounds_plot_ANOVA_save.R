peel_compounds_plot_save = function(comp.name,
                               file_for_metabo = "P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/peel/peel_files/files_for_metaboanalyst_18072022/pathway.name.for.metabo.18072022_for_ANOVA.csv",
                               path_for_save_file = "P:/yonatan/saved.pics.R/"
                               )
{
library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(stringr)
library(multcompView)
library(rlang)
library(ggpubr)

  # remove "/" if its in the end of "path_for_save_file" row
  if (substr(path_for_save_file, nchar(path_for_save_file)-1+1, nchar(path_for_save_file)) == "/") {
    path_for_save_file = substr(path_for_save_file,1,nchar(path_for_save_file)-1)
  }
  
# if there is a pthway name in the file from lc_pipline
  if (file.exists("P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/peel/peel_files/files_for_metaboanalyst_18072022/comp.and.pathway.table.18072022.csv") == T) {
    comp.and.pathway.file = read.csv("P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/peel/peel_files/files_for_metaboanalyst_18072022/comp.and.pathway.table.18072022.csv")
    pathway.name = comp.and.pathway.file[comp.and.pathway.file == comp.name,2]
    pathway.name = substring(pathway.name, 1,nchar(pathway.name)-3)
  } else {pathway.name = "plots"}
    
table = read.csv(file_for_metabo)
names(table) = table[1,]
names(table)[1] = "X"
table$X = gsub(" ", "_", table$X)

###
###

comp.name = comp.name
comp.name = gsub(" ", "_", comp.name)

###
###

title_x = ""  # or # names(data[1])
title_y = "Normalized Peak Area"  # or # names(data[2])
#title_y = ""

y.axis.hight = 1.1
y.label.hight = 1

bar.width = 0.5
error.width = 0.1

##############



data = data.frame(colnames(table), unlist(table[table$X == comp.name,], use.names = F))[-1,]
names(data) = data.frame(colnames(table), unlist(table[table$X == comp.name,], use.names = F))[1,]
data$X = gsub(" ","", data$X)
data$X = gsub("_",".", data$X)


data[,2] = as.numeric(data[,2])

title_main = names(data)[2]


X = names(data)[1]
Y = names(data[2])
level.order = unique(data$X)[c(6,1,7:10,2:5)]

###############
t.test_p.value = t.test(data[1:15,2],data[16:30,2], var.equal = T)$p.value
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

if (nchar(comp.name) >= 40) {
  comp.name = paste(substring(comp.name, 1,35),"....", substring(comp.name, nchar(comp.name)-1,nchar(comp.name)), sep = "")
  title_main = comp.name
}

#########################
stat = aov(as.formula(paste( Y,'~',X)), data) 
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

if (p_value_lable == "ns") {
  p_ns = annotate("text", x = 5.5, y = max(data[,2])*1.05, label= p_value_lable, size = 6)
} else {p_ns = annotate("text", x = 5.5, y = max(data[,2])*1.05, label= p_value_lable, size = 8)}

#########################
comp.plot = data %>%
  ggplot(aes(x=factor(x = !!sym(X), level=level.order), y = !!sym(Y))) +
#  geom_boxplot(fill=c(rep("grey70",5),rep("#562842",5))) + 
  geom_bar(stat = "summary",color="black" , fun.y = "mean", fill=c("#562842","grey70",rep("#562842",4), rep("grey70",4)), width=bar.width) +
  geom_errorbar (stat = "summary", fun.data = "mean_se", width = error.width) +

  theme_classic() +
  theme(axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 10, face="bold"),
        plot.title = element_text(size=10.35, face="bold.italic"),
        axis.text.y = element_text(size = 10, face="bold"),
        axis.line.x = element_line(size = 1.1),
        axis.line.y = element_line(size = 1.1),
        axis.ticks = element_line(size = 1.1) , 
        axis.ticks.length = unit(0.01, "cm")) +
  
  labs(title = title_main, x = title_x, y = title_y) + 
  
#   geom_jitter(color="#383838", size=0.625, alpha=0.9, width = 0.18)  +
#  stat_summary(fun.y=mean, geom="point", shape=23, size=4)  
  
  p_ns + 

 geom_text(data=LABELS, aes(x= X, y=(y.labels.order*y.label.hight), label=Letters), vjust = 0, fontface = "bold") + 

  scale_y_continuous(expand = c(0,0), limits = c(0,max(data[,2])*1.1))# +  ### change y.axis hight -> limits=c(0,hight))


dir.create(paste(path_for_save_file,"/lcms_peel_", Sys.Date(), sep = ""))

  dir.create(paste(path_for_save_file,"/lcms_peel_", Sys.Date(),"/",pathway.name, sep = ""))
  ggsave(paste(path_for_save_file,"/lcms_peel_", Sys.Date(),"/",pathway.name,"/", comp.name, ".png", sep = ""),
         plot = comp.plot, width = 3.64, height = 4.45)


print(title_main)
print(paste("p.value = ", round(t.test_p.value, digits = 15), sep = ""))
comp.plot
}
