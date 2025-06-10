# path and file name (csv)
csv_1 = "P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/cor.map/cor.map.2nd/comp.list.for.cor.map.updated.csv" # list of uniqe compounds and there pathways
csv_2 = "P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/LCMS_rnaseq_samples/LCMS_rnaseq.samples_190122.csv" # samples in order of lcms running
csv_3 = "P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/All_Compounds_pom_P_060722.csv" # data from lcms. need to be order like tmplate!





path.for.saving.files = "P:/yonatan/RNAseq_yonatan_2021/lcms/csv_files_for_metabo/comp.with.pathway.P/" # dont forget "/"  in the end of that row
date = "100722"





column_number = c(3,22:51)
first_row_fenotype = c("PUNI.marker", rep("high",15), rep("low",15))



to.save.files = T





#################################################

df.1 = read.csv(csv_1)[,1:2]
df.1 = df.1[!df.1$comp.pathway == "",]

#samples = read.csv("../LCMS_rnaseq.samples_190122.csv")
samples = read.csv(csv_2)

lc.data = read.csv(csv_3)
lc.mzcl = lc.data[!lc.data$ï..Tags == "",]
lc = lc.mzcl[,column_number] # transgenic
lc.for.met = lc
names(lc.for.met) = c("Name", samples$sample[as.numeric(stringr::str_extract(names(lc.for.met), "\\d+")[-1])])
lc.for.met$Name = gsub("std", "STD", lc.for.met$Name)

# merge to create table for.metabo and for names.with.pathway
comp.pathway.name.table.with.metabo = merge.data.frame(lc.for.met,df.1, by = "Name", all.y = T)

# add puni.marker to the 1st row
comp.pathway.name.table.with.metabo[nrow(comp.pathway.name.table.with.metabo)+1,1:31] = first_row_fenotype
comp.pathway.name.table.with.metabo = comp.pathway.name.table.with.metabo[c(nrow(comp.pathway.name.table.with.metabo),
                                                                            1:nrow(comp.pathway.name.table.with.metabo)-1),]

# make non-duplicates "Name"
string1 = comp.pathway.name.table.with.metabo$Name
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
comp.pathway.name.table.with.metabo$Name = mstring1 

# make non-duplicates "comp.pathway"
string2 = comp.pathway.name.table.with.metabo$comp.pathway
mstring2 <- make.unique(as.character(string2), sep=" X" )
tmp <- !duplicated(string2)
for (i in 1:length(mstring2[tmp])){
  mstring2[tmp][i]<-ifelse(string2[tmp][i] %in% string2[duplicated(string2)]
                           , gsub("(.*)","\\1 X0", mstring2[tmp][i])
                           , mstring2[tmp][i]
  )
}
end <- sub(".* X([0-9]+)","\\1",grep(" X([0-9]*)$",mstring2,value=T) )
beg <- sub("(.* X)[0-9]+","\\1",grep(" X([0-9]*)$",mstring2,value=T) )
newend <- as.numeric(end)+1
mstring2[grep(" X([0-9]*)$",mstring2)] <- paste0(beg,newend)
comp.pathway.name.table.with.metabo$comp.pathway = mstring2

# save files
file1 = comp.pathway.name.table.with.metabo[,c(length(comp.pathway.name.table.with.metabo),2:(length(comp.pathway.name.table.with.metabo)-1))]
file1[1,1] = comp.pathway.name.table.with.metabo[1,1]

if (to.save.files == T) {
  write.csv(file1,
            paste(path.for.saving.files,"pathway.name.for.metabo.",date,".csv", sep = ""),
            row.names = F)
  write.csv(comp.pathway.name.table.with.metabo[,1:(length(comp.pathway.name.table.with.metabo)-1)],
            paste(path.for.saving.files,"comp.name.for.metabo.",date,".csv", sep = ""),
            row.names = F)
  write.csv(comp.pathway.name.table.with.metabo[-1,c(1,length(comp.pathway.name.table.with.metabo))],
            paste(path.for.saving.files,"comp.and.pathway.table.",date,".csv", sep = ""),
            row.names = F)
}

