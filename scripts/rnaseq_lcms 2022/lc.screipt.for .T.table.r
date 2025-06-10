setwd("P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/")

samples = read.csv("../LCMS_rnaseq.samples_190122.csv")

lc.data = read.csv("All_Compounds_pom_T_060722.csv")
lc.mzcl = lc.data[!lc.data$ï..Tags == "",]
#lc = lc.mzcl[,c(4,17:46)] # peel
lc = lc.mzcl[,c(3,25:44)] # transgenic

lc.for.met = lc
lc.for.met[nrow(lc)+1,] = c("OE","bHLH94","bHLH94","bHLH94","bHLH94","bHLH94",
                            "EV","EV","EV","EV","EV",
                            "MYB6","MYB6","MYB6","MYB6","MYB6",
                            "MYB8b","MYB8b","MYB8b","MYB8b","MYB8b")
lc.for.met = lc.for.met[c(nrow(lc)+1,1:nrow(lc)),]
names(lc.for.met)[c(7,16)] = c("Norm..Area..4.raw..F2.","Norm..Area..14.raw..F12.")
#lc.for.met[1,1] = "OE"
#lc.for.met[1,17:31] = "low"

names(lc.for.met) = c("Name", samples$sample[as.numeric(stringr::str_extract(names(lc.for.met), "\\d+")[-1])])
#lc.for.met = lc.for.met[!duplicated(lc.for.met$Name),]

lc.for.met.save = lc.for.met[,c(1,7:11,2:6,12:21)]

lc.for.met.save$Name = gsub("std", "STD", lc.for.met.save$Name)

string = lc.for.met.save$Name
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
lc.for.met.save$Name = mstring 

write.csv(lc.for.met.save,"../csv_files_for_metabo/lc.T.for.metaboanalyst.csv", row.names = F)
