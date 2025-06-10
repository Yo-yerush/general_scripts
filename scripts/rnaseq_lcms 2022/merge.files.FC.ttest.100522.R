setwd("C:/Users/yonatany/Desktop/100522/merge files")
FC = read.csv("fold_change bhlh_VS_ev.csv")
ttest = read.csv("t_test bhlh_VS_ev.csv")

new.table = merge.data.frame(FC, ttest, by = "X", all.x = T)

new.table2 = new.table[,c(1,3:5)]

new.table2 = new.table2[order(new.table2$log2.FC.,decreasing = T),]

write.csv(new.table2,"FC_ttest_table_bhlh.VS.ev_100522.csv" ,row.names = F)
