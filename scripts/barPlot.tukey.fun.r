barPlot.tukey = function(data, title_x = names(data[1]), title_y = names(data[2]), title_main = "", y.axis.hight = max(y.labels.order)*1.1, label.hight = y.labels.order, bar.width = 0.5, error.width = 0.1, jitter = F)
{
  library(dplyr)
  library(ggplot2)
  library(ggsignif)
  library(tidyr)
  library(stringr)
  library(multcompView)
  library(rlang)
  
  ###
  X = names(data[1])
  Y = names(data[2])
  level.order = unique(data[,1])
  
  ##########################
  
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
  
  if (jitter == T) {
    p.jitter = geom_jitter(color="steelblue4", size=1, alpha=0.9, width = 0.18)
  } else {p.jitter = geom_blank()}
  
  data %>%
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
    scale_y_continuous(expand = c(0,0), limits = c(0,y.axis.hight)) +  ### change y.axis hight -> limits=c(0,hight))
    geom_text(data=LABELS, aes(x= X, y=(label.hight), label=Letters), vjust = 0, fontface = "bold") +  ### change tukey labels hight - y=(hight)
    p.jitter
}