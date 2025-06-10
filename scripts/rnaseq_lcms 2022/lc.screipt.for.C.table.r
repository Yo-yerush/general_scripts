
samples = read.csv("P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/LCMS_rnaseq.samples_C_260522.csv")

lc.data = read.csv("P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/first.table/Compounds_C_290522.csv")
#lc.mzcl = lc.data[!lc.data$ï..Tags == "",]

lc = lc.data[,c(3,23,25,29:31,26,27,35:37)] # clone

lc.for.met = lc
lc.for.met[nrow(lc)+1,] = c("PUNI.marker", rep("high",2), rep("low",3), rep("high",2), rep("low",3))
lc.for.met = lc.for.met[c(nrow(lc)+1,1:nrow(lc)),]

#lc.for.met[1,1] = "OE"
#lc.for.met[1,17:31] = "low"

names(lc.for.met) = c("Name", rep("101-1",2), rep("116-17",3), rep("Ban-G",2), rep("Dan-A",3))

lc.for.met$Name = gsub("std", "STD", lc.for.met$Name)




# add X to duplicated rows
string = lc.for.met$Name
mstring <- make.unique(as.character(string), sep=" X" )
tmp <- !duplicated(string)
for (i in 1:length(mstring[tmp])){
  mstring[tmp][i]<-ifelse(string[tmp][i] %in% string[duplicated(string)]
                          , gsub("(.*)","\\1 X0", mstring[tmp][i])
                          , mstring[tmp][i]
  )
}
end <- sub(".* X([0-9]+)","\\1",grep(" X([0-9]*)$",mstring,value=T) )
beg <- sub("(.* X)[0-9]+","\\1",grep(" X([0-9]*)$",mstring,value=T) )
newend <- as.numeric(end)+1
mstring[grep(" X([0-9]*)$",mstring)] <- paste0(beg,newend)
lc.for.met$Name = mstring 



write.csv(lc.for.met,"P:/pomegranate RNAseq 2022/clone/lcms.files/pom.C.for.metaboanalyst.csv", row.names = F)



