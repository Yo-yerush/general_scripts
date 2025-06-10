comp_pathway = read.csv("P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/transgenic/transgenic_files/files_for_metaboanalyst_18072022/comp.and.pathway.table.18072022.csv")
fc = read.csv("P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/transgenic/metabo_outputs/fold_change_all.csv")
ttest = read.csv("P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/transgenic/metabo_outputs/t_test_all.csv")[,1:3]
names(fc)[1] = names(comp_pathway)[1]
names(ttest)[1] = names(comp_pathway)[1]

fc_t = merge.data.frame(ttest, fc, by = names(fc)[1])
all = merge.data.frame(fc_t, comp_pathway, by = names(fc)[1])

all$comp.pathway = gsub(" X[0-9]", "", all$comp.pathway)
all$comp.pathway = gsub("[0-9]", "", all$comp.pathway)

#write.csv(all, "P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/transgenic/metabo_outputs/t_test_FC_with_pathway_2.csv", row.names = F)

lc.data = read.csv("P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/All_Compounds_pom_T_060722.csv")
names(lc.data)[1] = "Tags"
lc.data = lc.data[!lc.data$Tags == "",]
lc.data$Name = gsub("std", "STD", lc.data$Name)
lc.data = lc.data[,c(3,4,11,13,14,23)]
names(lc.data)[3:5] = c("Molecular Weight","M/Z", "RT")

string1 = lc.data$Name
mstring1 <- make.unique(as.character(string1), sep=" X" )
tmp <- !duplicated(string1)
for (i in 1:length(mstring1[tmp])){
  mstring1[tmp][i]<-ifelse(string1[tmp][i] %in% string1[duplicated(string1)]
                           , gsub("(.*)","\\1 X0", mstring1[tmp][i])
                           , mstring1[tmp][i]
  )
}
end <- sub(".* X([0-9]+)","\\1",grep(" X([0-9]*)$",mstring1,value=T) )
beg <- sub("(.* X)[0-9]+","\\1",grep(" X([0-9]*)$",mstring1,value=T) )
newend <- as.numeric(end)+1
mstring1[grep(" X([0-9]*)$",mstring1)] <- paste0(beg,newend)
lc.data$Name = mstring1 


lcms_table = merge.data.frame(all, lc.data, by = names(all)[1])
names(lcms_table)[c(3:6,11)] = c("p.value","Fold Change", "log2(Fold Change)", "Class", "MS2")
lcms_table$MS2 = gsub("1","",lcms_table$MS2)
lcms_table$MS2 = gsub("2","*",lcms_table$MS2)

write.csv(lcms_table, "P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/transgenic/lcms_table_bHLH94.csv", row.names = F)
