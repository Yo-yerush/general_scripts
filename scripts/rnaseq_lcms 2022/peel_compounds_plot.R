peel_compounds_plot = function(comp.name, 
                               file_for_metabo = "P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/peel/peel_files/files_for_metaboanalyst_18072022/comp.name.for.metabo.18072022.csv")
{
library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(stringr)
library(multcompView)
library(rlang)
library(ggpubr)


table = read.csv(file_for_metabo)
names(table) = table[1,]
names(table)[1] = "X"



###
###

comp.name = comp.name


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

########
#names(data) = c("data1","data2")
#df_p_val <- rstatix::t_test(data, data2 ~ data1, ref.group = "EV") %>% 
#  rstatix::add_xy_position()
########

names(data) = data.frame(colnames(table), unlist(table[table$X == comp.name,], use.names = F))[1,]
data[,2] = as.numeric(data[,2])

title_main = names(data)[2]


X = names(data)[1]
Y = names(data[2])
level.order = unique(data$X)

###############
t.test_p.value = t.test(data[1:15,2],data[16:30,2], var.equal = T)$p.value
if (t.test_p.value <= 0.05 & t.test_p.value > 0.01) {
  p_value_lable = "*"
} else if (t.test_p.value <= 0.01 & t.test_p.value > 0.001) {
  p_value_lable = "**"
} else if (t.test_p.value <= 0.001) {
  p_value_lable = "***"
}

#########################


#c("firebrick4","grey70")
#########################
comp.plot = data %>%
  ggplot(aes(x=factor(x = !!sym(X), level=level.order), y = !!sym(Y))) +
  geom_boxplot(fill=c("grey70","#562842"), ) + 

  theme_classic() +
  theme(axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 10, face="bold"),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.text.y = element_text(size = 10, face="bold"),
        axis.line.x = element_line(size = 1.1),
        axis.line.y = element_line(size = 1.1),
        axis.ticks = element_line(size = 1.1) , 
        axis.ticks.length = unit(0.01, "cm")) +
  
  labs(title = title_main, x = title_x, y = title_y) + 
  
   geom_jitter(color="#383838", size=2, alpha=0.9, width = 0.18)  +
#  stat_summary(fun.y=mean, geom="point", shape=23, size=4)  
  
  annotate("text", x = 1.5, y = max(data[,2])*1, label= p_value_lable, size = 7)


#  scale_y_continuous(expand = c(0,0), limits = c(0,max(data[,2])*1.1))# +  ### change y.axis hight -> limits=c(0,hight))
print(title_main)
print(paste("p.value = ", round(t.test_p.value, digits = 15), sep = ""))
comp.plot
}
