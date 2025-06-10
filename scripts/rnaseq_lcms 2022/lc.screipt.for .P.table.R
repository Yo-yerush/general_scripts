setwd("P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/")

samples = read.csv("../LCMS_rnaseq.samples_190122.csv")

lc.data = read.csv("All_Compounds_pom_P_060722.csv")
lc.mzcl = lc.data[!lc.data$ï..Tags == "",]
#lc = lc.mzcl[,c(4,17:46)] # peel
lc = lc.mzcl[,c(3,22:51)] # transgenic

lc.for.met = lc
lc.for.met[nrow(lc)+1,] = c("PUNI.marker", rep("high",15), rep("low",15))
lc.for.met = lc.for.met[c(nrow(lc)+1,1:nrow(lc)),]

#lc.for.met[1,1] = "OE"
#lc.for.met[1,17:31] = "low"

names(lc.for.met) = c("Name", samples$sample[as.numeric(stringr::str_extract(names(lc.for.met), "\\d+")[-1])])
#lc.for.met = lc.for.met[!duplicated(lc.for.met$Name),]


lc.for.met$Name = gsub("std", "STD", lc.for.met$Name)

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

write.csv(lc.for.met,"../csv_files_for_metabo/lc.P.for.metaboanalyst.csv", row.names = F)
