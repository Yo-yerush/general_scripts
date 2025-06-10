# Ctrl+Alt+T runs the current code section
# library ----------------------------------------------------------
# install.packages("ggplot2")
# install.packages("cowplot")
# install.packages("multcompView")
# install.packages("multcomp")
# install.packages("agricolae")
# install.packages("stringr")


# load data --------------------------------------------------------

ds1 = read.csv(choose.files()) # has to be in CSV format

########################################################################
# User input !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    #
                                                                      #
# column number input is needed !!!                                   #
gene_name = as.integer(1) # Enter 'gene name' column number           #
padj = as.integer(6) # Enter 'padj' column number                     #
GO_name = as.integer(26) # Enter 'GO name' column number              #
ctrl_samp = (7:9)
treatment_samp = c(10,11,13,14)
cut_off_ = 10 # the top 'X' changed GO's                              #
saving_pathway = "M:\\1.Lior Lab\\Shay\\RNASEQ\\R_output"
padj_val = as.numeric(0.05) # Enter 'padj' cutoff value
#######################################################################

library(ggplot2)
library(cowplot)
library(multcompView)
library(multcomp)
library(agricolae)
library(stringr)
library(writexl)


ds1 = ds1[which(ds1$padj <= padj_val),] # subsisting 

ds1.2 = ds1[,c(gene_name,padj,GO_name,ctrl_samp,treatment_samp)]
colnames(ds1.2)[1:3] = c("gene_name","padj","GO_name")
ds1.2 = ds1.2[which(ds1.2$padj <= padj_val),] # subsisting 

ds1.2$ctrl_mean = apply(ds1[,ctrl_samp], 1, mean)
ds1.2$s_mean = apply(ds1[,treatment_samp], 1, mean)

number_changed = dim(ds1.2)[1]
number_up = sum(ds1.2$ctrl_mean < ds1.2$s_mean)
number_down = sum(ds1.2$ctrl_mean > ds1.2$s_mean)

ds_gene_num = as.data.frame(c(number_changed,number_up,number_down))
colnames(ds_gene_num) = "freq"
row.names(ds_gene_num)

label_angel = 0
margins = unit(c(0.5,0.5,0.5,0.5),"cm")# plot.margin = unit(c(TOP,RIGHT,BOTTOM,LEFT),"cm")
y_lim = c(0,2000)


plot_gene_num = ggplot(data = ds_gene_num, aes(x = c("Changed (up or down)","Upregulated","Down regulated"), y = freq)) +
  geom_bar(stat = "identity") + geom_text(data = ds_gene_num, aes(label = freq),position = position_nudge(y = 50)) +
  theme(axis.text.x = element_text(angle = label_angel,hjust = 0.5, vjust = 0.5),
        legend.position = "none", plot.margin = margins) + ggtitle("segnificantly different expressed genes") + 
  ylim(y_lim) + theme(plot.title = element_text(hjust = 0.5)) + labs(y="Number of genes", x="")
plot_gene_num

# generation down and up regulated files ----------------------------------------------------

ds2_s_down = ds1.2[which(ds1.2$ctrl_mean > ds1.2$s_mean),]
ds2_s_up = ds1.2[which(ds1.2$ctrl_mean < ds1.2$s_mean),]

#  ALL (changed) DS --------------------------------------

GO_names = strsplit(x =  ds1.2$GO_name,split =  "; ") #for splitting GO names
GO_names = unlist(GO_names)

GO_names2 = substr(x = GO_names, start = 0, stop = 1) #
GO_names2 = factor(GO_names2)
GO_names3 = substr(x = GO_names, start = 3, stop = 50) #
GO_names3 = factor(GO_names3)
cbind(GO_names2,GO_names3)

ds2 = data.frame("GO_cat" = GO_names2, "GO_nam" = GO_names3)

go_categories_names = c("Cellular Component", "Molecular Function", "Biological Process")

ds2$GO_cate2 = NA

ds2$GO_cate2[which(ds2$GO_cat == "C")] = "Cellular Component"
ds2$GO_cate2[which(ds2$GO_cat == "F")] = "Molecular Function"
ds2$GO_cate2[which(ds2$GO_cat == "P")] = "Biological Process"

ds2$GO_cate2 = factor(ds2$GO_cate2)

ds2_s_changed = ds2

ds2 = NULL
#  up regulated DS --------------------------------------

GO_names = strsplit(x =  ds2_s_up$GO_name,split =  "; ") #for splitting GO names
GO_names = unlist(GO_names)

GO_names2 = substr(x = GO_names, start = 0, stop = 1) #
GO_names2 = factor(GO_names2)
GO_names3 = substr(x = GO_names, start = 3, stop = 50) #
GO_names3 = factor(GO_names3)
cbind(GO_names2,GO_names3)

ds2 = data.frame("GO_cat" = GO_names2, "GO_nam" = GO_names3)

go_categories_names = c("Cellular Component", "Molecular Function", "Biological Process")

ds2$GO_cate2 = NA

ds2$GO_cate2[which(ds2$GO_cat == "C")] = "Cellular Component"
ds2$GO_cate2[which(ds2$GO_cat == "F")] = "Molecular Function"
ds2$GO_cate2[which(ds2$GO_cat == "P")] = "Biological Process"

ds2$GO_cate2 = factor(ds2$GO_cate2)

ds2_s_up = ds2

ds2 = NULL
#  down regulated DS --------------------------------------

GO_names = strsplit(x =  ds2_s_down$GO_name,split =  "; ") #for splitting GO names
GO_names = unlist(GO_names)

GO_names2 = substr(x = GO_names, start = 0, stop = 1) #
GO_names2 = factor(GO_names2)
GO_names3 = substr(x = GO_names, start = 3, stop = 50) #
GO_names3 = factor(GO_names3)
cbind(GO_names2,GO_names3)

ds2 = data.frame("GO_cat" = GO_names2, "GO_nam" = GO_names3)

go_categories_names = c("Cellular Component", "Molecular Function", "Biological Process")

ds2$GO_cate2 = NA

ds2$GO_cate2[which(ds2$GO_cat == "C")] = "Cellular Component"
ds2$GO_cate2[which(ds2$GO_cat == "F")] = "Molecular Function"
ds2$GO_cate2[which(ds2$GO_cat == "P")] = "Biological Process"

ds2$GO_cate2 = factor(ds2$GO_cate2)

ds2_s_down = ds2

# CHANGED  ###############################################################
C_C = subset(ds2_s_changed, subset = GO_cate2 == "Cellular Component")
C_C = table(C_C$GO_nam)
C_C =  C_C[C_C != 0]
C_C = C_C[order(C_C, decreasing = TRUE)]

M_F = subset(ds2_s_changed, subset = GO_cate2 == "Molecular Function")
M_F = table(M_F$GO_nam)
M_F =  M_F[M_F != 0]
M_F = M_F[order(M_F, decreasing = TRUE)]

B_P = subset(ds2_s_changed, subset = GO_cate2 == "Biological Process")
B_P = table(B_P$GO_nam)
B_P =  B_P[B_P != 0]
B_P = B_P[order(B_P, decreasing = TRUE)]

GO_freq = c(C_C,M_F,B_P)
# subset only high freq ------------------------------
cut_off = c(1:cut_off_) #which(C_C > 150)

C_C2_changed = C_C[cut_off]
C_C2_changed = as.data.frame(C_C2_changed)
colnames(C_C2_changed) = c("GO_NAME","Freq")

M_F2_changed = M_F[cut_off]
M_F2_changed = as.data.frame(M_F2_changed)
colnames(M_F2_changed) = c("GO_NAME","Freq")

B_P2_changed = B_P[cut_off]
B_P2_changed = as.data.frame(B_P2_changed)
colnames(B_P2_changed) = c("GO_NAME","Freq")

max_changed = c(max(C_C),max(M_F),max(B_P))
# UP  ###############################################################
C_C = subset(ds2_s_up, subset = GO_cate2 == "Cellular Component")
C_C = table(C_C$GO_nam)
C_C =  C_C[C_C != 0]
C_C = C_C[order(C_C, decreasing = TRUE)]

M_F = subset(ds2_s_up, subset = GO_cate2 == "Molecular Function")
M_F = table(M_F$GO_nam)
M_F =  M_F[M_F != 0]
M_F = M_F[order(M_F, decreasing = TRUE)]

B_P = subset(ds2_s_up, subset = GO_cate2 == "Biological Process")
B_P = table(B_P$GO_nam)
B_P =  B_P[B_P != 0]
B_P = B_P[order(B_P, decreasing = TRUE)]

GO_freq = c(C_C,M_F,B_P)
# subset only high freq ------------------------------
cut_off = c(1:cut_off_) #which(C_C > 150)

C_C2_up = C_C[cut_off]
C_C2_up = as.data.frame(C_C2_up)
colnames(C_C2_up) = c("GO_NAME","Freq")

M_F2_up = M_F[cut_off]
M_F2_up = as.data.frame(M_F2_up)
colnames(M_F2_up) = c("GO_NAME","Freq")

B_P2_up = B_P[cut_off]
B_P2_up = as.data.frame(B_P2_up)
colnames(B_P2_up) = c("GO_NAME","Freq")

max_up = c(max(C_C),max(M_F),max(B_P))
# DOWN  ###############################################################
C_C = subset(ds2_s_down, subset = GO_cate2 == "Cellular Component")
C_C = table(C_C$GO_nam)
C_C =  C_C[C_C != 0]
C_C = C_C[order(C_C, decreasing = TRUE)]

M_F = subset(ds2_s_down, subset = GO_cate2 == "Molecular Function")
M_F = table(M_F$GO_nam)
M_F =  M_F[M_F != 0]
M_F = M_F[order(M_F, decreasing = TRUE)]

B_P = subset(ds2_s_down, subset = GO_cate2 == "Biological Process")
B_P = table(B_P$GO_nam)
B_P =  B_P[B_P != 0]
B_P = B_P[order(B_P, decreasing = TRUE)]

GO_freq = c(C_C,M_F,B_P)
# subset only high freq ------------------------------
cut_off = c(1:cut_off_) #which(C_C > 150)

C_C2_down = C_C[cut_off]
C_C2_down = as.data.frame(C_C2_down)
colnames(C_C2_down) = c("GO_NAME","Freq")

M_F2_down = M_F[cut_off]
M_F2_down = as.data.frame(M_F2_down)
colnames(M_F2_down) = c("GO_NAME","Freq")

B_P2_down = B_P[cut_off]
B_P2_down = as.data.frame(B_P2_down)
colnames(B_P2_down) = c("GO_NAME","Freq")

max_down = c(max(C_C),max(M_F),max(B_P))
# plotting -----------------------------------------

# changed  ###############################################################
# plot 1 ---------------------------------------

label_angel = 290
margins = unit(c(0.5,3,0.5,0.5),"cm")# plot.margin = unit(c(TOP,RIGHT,BOTTOM,LEFT),"cm")
lim_y = c(0,round((max_changed[1]/100)+0.5)*100)

plot1 = ggplot(data = C_C2_changed, aes(x = reorder(GO_NAME, -Freq) ,y = Freq)) +
  geom_bar( stat = "identity") + 
  theme(axis.text.x = element_text( angle = label_angel,hjust = 0, vjust = 0.5),
        legend.position = "none", plot.margin = margins) + labs(x = "") +
  ggtitle("Cellular Component - changed") + ylim(lim_y) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Freq),position = position_nudge(y = (lim_y[2]/40)))

# plot 2 ---------------------------------------

label_angel = 290
margins = unit(c(0.5,3,0.5,0.5),"cm")# plot.margin = unit(c(TOP,RIGHT,BOTTOM,LEFT),"cm")
lim_y = c(0,round((max_changed[2]/100)+0.5)*100)

plot2 = ggplot(data = M_F2_changed, aes(x = reorder(GO_NAME, -Freq) ,y = Freq)) +
  geom_bar( stat = "identity") + 
  theme(axis.text.x = element_text( angle = label_angel,hjust = 0, vjust = 0.5),
        legend.position = "none", plot.margin = margins,
        axis.title.y = element_blank()) + labs(x = "")+
  ggtitle("Molecular Function - changed") + ylim(lim_y) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Freq),position = position_nudge(y = (lim_y[2]/40)))

# plot 3 ---------------------------------------

label_angel = 290
margins = unit(c(0.5,3,0.5,0.5),"cm")# plot.margin = unit(c(TOP,RIGHT,BOTTOM,LEFT),"cm")
lim_y = c(0,round((max_changed[3]/100)+0.5)*100)

plot3 = ggplot(data = B_P2_changed, aes(x = reorder(GO_NAME, -Freq) ,y = Freq)) +
  geom_bar( stat = "identity") + 
  theme(axis.text.x = element_text( angle = label_angel,hjust = 0, vjust = 0.5),
        legend.position = "none", plot.margin = margins,
        axis.title.y = element_blank()) + labs(x = "")+
  ggtitle("Biological Process - changed")  + ylim(lim_y) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Freq),position = position_nudge(y = (lim_y[2]/40)))


# down ###############################################################
# plot 1 ---------------------------------------

label_angel = 290
margins = unit(c(0.5,3,0.5,0.5),"cm")# plot.margin = unit(c(TOP,RIGHT,BOTTOM,LEFT),"cm")
lim_y = c(0,round((max_down[1]/100)+0.5)*100)

plot4 = ggplot(data = C_C2_down, aes(x = reorder(GO_NAME, -Freq) ,y = Freq)) +
geom_bar( stat = "identity") + 
  theme(axis.text.x = element_text( angle = label_angel,hjust = 0, vjust = 0.5),
        legend.position = "none", plot.margin = margins) + labs(x = "") +
  ggtitle("Cellular Component - down") + ylim(lim_y) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Freq),position = position_nudge(y = (lim_y[2]/40)))

# plot 2 ---------------------------------------

label_angel = 290
margins = unit(c(0.5,3,0.5,0.5),"cm")# plot.margin = unit(c(TOP,RIGHT,BOTTOM,LEFT),"cm")
lim_y = c(0,round((max_down[2]/100)+0.5)*100)

plot5 = ggplot(data = M_F2_down, aes(x = reorder(GO_NAME, -Freq) ,y = Freq)) +
  geom_bar( stat = "identity") + 
  theme(axis.text.x = element_text( angle = label_angel,hjust = 0, vjust = 0.5),
        legend.position = "none", plot.margin = margins,
        axis.title.y = element_blank()) + labs(x = "")+
  ggtitle("Molecular Function - down") + ylim(lim_y) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Freq),position = position_nudge(y = (lim_y[2]/40)))

# plot 3 ---------------------------------------

label_angel = 290
margins = unit(c(0.5,3,0.5,0.5),"cm")# plot.margin = unit(c(TOP,RIGHT,BOTTOM,LEFT),"cm")
lim_y = c(0,round((max_down[3]/100)+0.5)*100)

plot6 = ggplot(data = B_P2_down, aes(x = reorder(GO_NAME, -Freq) ,y = Freq)) +
  geom_bar( stat = "identity") + 
  theme(axis.text.x = element_text( angle = label_angel,hjust = 0, vjust = 0.5),
        legend.position = "none", plot.margin = margins,
        axis.title.y = element_blank()) + labs(x = "")+
  ggtitle("Biological Process - down")  + ylim(lim_y) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Freq),position = position_nudge(y = (lim_y[2]/40)))

# up  ###############################################################
# plot 1 ---------------------------------------

label_angel = 290
margins = unit(c(0.5,3,0.5,0.5),"cm")# plot.margin = unit(c(TOP,RIGHT,BOTTOM,LEFT),"cm")
y_lim = c(0,round((max_up[1]/100)+0.5)*100)

plot7 = ggplot(data = C_C2_up, aes(x = reorder(GO_NAME, -Freq) ,y = Freq)) +
  geom_bar( stat = "identity") + 
  theme(axis.text.x = element_text( angle = label_angel,hjust = 0, vjust = 0.5),
        legend.position = "none", plot.margin = margins) + labs(x = "") +
  ggtitle("Cellular Component - up") + ylim(y_lim) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Freq),position = position_nudge(y = (y_lim[2]/40)))

# plot 2 ---------------------------------------

label_angel = 290
margins = unit(c(0.5,3,0.5,0.5),"cm")# plot.margin = unit(c(TOP,RIGHT,BOTTOM,LEFT),"cm")
y_lim = c(0,round((max_up[2]/100)+0.5)*100)


plot8 = ggplot(data = M_F2_up, aes(x = reorder(GO_NAME, -Freq) ,y = Freq)) +
  geom_bar( stat = "identity") + 
  theme(axis.text.x = element_text( angle = label_angel,hjust = 0, vjust = 0.5),
        legend.position = "none", plot.margin = margins,
        axis.title.y = element_blank()) + labs(x = "")+
  ggtitle("Molecular Function - up") + ylim(y_lim) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Freq),position = position_nudge(y = (y_lim[2]/40)))

# plot 3 ---------------------------------------

label_angel = 290
margins = unit(c(0.5,3,0.5,0.5),"cm")# plot.margin = unit(c(TOP,RIGHT,BOTTOM,LEFT),"cm")
y_lim = c(0,round((max_up[3]/100)+0.5)*100)


plot9 = ggplot(data = B_P2_up, aes(x = reorder(GO_NAME, -Freq) ,y = Freq)) +
  geom_bar( stat = "identity") + 
  theme(axis.text.x = element_text( angle = label_angel,hjust = 0, vjust = 0.5),
        legend.position = "none", plot.margin = margins,
        axis.title.y = element_blank()) + labs(x = "")+
    ggtitle("Biological Process - up") + ylim(y_lim) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Freq),position = position_nudge(y = (y_lim[2]/40)))

#--------------------------------------------------------------------------

plot_grid(plot1,plot2,plot3, align = "hv",  nrow = 1 , rel_widths =)
ggsave(
  "GO_changed.png",
  path = saving_pathway ,
  scale = 1,
  width = 30,
  height = 15,
  units = c("cm"),
  dpi = 300,
)
plot_grid(plot4,plot5,plot6,align = "hv",  nrow = 1 , rel_widths =)
ggsave(
  "GO_down.png",
  path = saving_pathway ,
  scale = 1.5,
  width = 30,
  height = 15,
  units = c("cm"),
  dpi = 300,
)
plot_grid(plot7,plot8,plot9,align = "hv",  nrow = 1 , rel_widths =)
ggsave(
  "GO_up.png",
  path = saving_pathway ,
  scale = 1.5,
  width = 30,
  height = 15,
  units = c("cm"),
  dpi = 300,
)

plot_gene_num
ggsave(
  "changed_genes_count.png",
  path = saving_pathway ,
  scale = 1.5,
  width = 10,
  height = 15,
  units = c("cm"),
  dpi = 300,
)
#--------------------------------------------
ds2_go_changed = table(ds2_s_changed$GO_nam)
ds2_go_changed = as.data.frame(ds2_go_changed)

ds2_go_changed = ds2_go_changed[order(ds2_go_changed$Freq, decreasing = TRUE),]
write_xlsx(ds2_go_changed,paste(saving_pathway,"\\changed.xlsx",sep = ""))
####################################################
ds2_go_up = table(ds2_s_up$GO_nam)
ds2_go_up = as.data.frame(ds2_go_up)

ds2_go_up = ds2_go_up[order(ds2_go_up$Freq, decreasing = TRUE),]
write_xlsx(ds2_go_up,paste(saving_pathway,"\\up.xlsx",sep = ""))
####################################################
ds2_go_down = table(ds2_s_down$GO_nam)
ds2_go_down = as.data.frame(ds2_go_down)

ds2_go_down = ds2_go_down[order(ds2_go_down$Freq, decreasing = TRUE),]
write_xlsx(ds2_go_down,paste(saving_pathway,"\\down.xlsx",sep = ""))

