prepare_metabo_table <- function(table.from.lcms.software,
                                  lcms.samples.file, path.for.saving.files, 
                                  group_title_row, sample_order_name = F,
                                  Tags = T) {


##############
  
### table.from.lcms.software
# path to raw data file from lcms (csv file from lcms software in Tel-Hai)
#   example: table.from.lcms.software = "P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/All_Compounds_pom_P_060722.csv"

  
### lcms.samples.file
# path to samples file (csv) in order of lcms running, with BLK and QC
#   example: lcms.samples.file = "P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/LCMS_rnaseq_samples/LCMS_rnaseq.samples_190122.csv"

#############

  
### path.for.saving.files 
# path for the saved files directory
#   example: path.for.saving.files = "P:/yonatan/RNAseq_yonatan_2021/lcms/csv_files_for_metabo/comp.with.pathway.P/" 



### group_title_row
# group description for metaboanalyst analysis
#   example for 30 samples and 3 groups: 
# group_title_row = c("OE_plants", rep("wild.type",10), rep("mt.line1",10), rep("mt.line2",10))

### sample_order_name 
# to get sample names easily, run script with: 
#   all  non-default arguments except <group_title_row> and with <sample_order_name = T>
#   *no file will saved

  
### Tags
# logical, if used "Tags" column for filtering in the lcms software

##############
  
  
  date = gsub("/","",format(Sys.time(), "%x"))
  
# remove "/" if its in the end of "path.for.saving.files" row
if (substr(path.for.saving.files, nchar(path.for.saving.files)-1+1, nchar(path.for.saving.files)) == "/") {
  path.for.saving.files = substr(path.for.saving.files,1,nchar(path.for.saving.files)-1)
}

# samples file
samples = read.csv(lcms.samples.file)

# lcms raw data file
lc.data = read.csv(table.from.lcms.software)
names(lc.data)[1] = "Tags"

if (Tags == T) {
  lc.mzcl = lc.data[!lc.data$Tags == "",]
} else {
  lc.mzcl = lc.data
}

lc = lc.mzcl[,c(which(names(lc.data) == "Name"),
                grep("Norm..Area..[0-9].raw|Norm..Area..[0-9][0-9].raw|Norm..Area..[0-9][0-9][0-9].raw|Norm..Area..[0-9]_[0-9].raw", 
                   names(lc.data), ignore.case = TRUE, value = F))]

for (i in 1:length(lc[,1])) {
  lc[i,] = gsub("\r|\n","",lc[i,])
}

lc.for.met = lc
names(lc.for.met) = c("Name", samples$sample[as.numeric(stringr::str_extract(names(lc.for.met),
                                                                             "\\d+")[-1])])
lc.for.met$Name = gsub("std", "STD", lc.for.met$Name)


# if you want to get only <sample_order_name>
if (sample_order_name == T) {
  print(names(lc.for.met)[-1])
     opt <- options(show.error.messages = FALSE)
     on.exit(options(opt))
     stop()
}


# add a group name to the 1st row
lc.for.met[(nrow(lc.for.met)+1),1:(length(lc.for.met))] = group_title_row
lc.for.met = lc.for.met[c(nrow(lc.for.met), 
                          1:(nrow(lc.for.met)-1)),]

# make non-duplicates "Name"
string1 = lc.for.met$Name
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
lc.for.met$Name = mstring1 


# save files
if (file.exists(paste(path.for.saving.files,"/table_for_metaboanalyst_",date,".csv", sep = "")) == T) {
  message(paste("the file 'table_for_metaboanalyst_",date,
                             "' is already exist, do you want to continue?\n  *it will remove the old file with the same name in that folder", sep = ""))
  
  if (regexpr(readline("\t(y/n): "), 'y', ignore.case = TRUE) == 1) {
    
    write.csv(lc.for.met,
              paste(path.for.saving.files,"/table_for_metaboanalyst_",date,".csv", sep = ""),
              row.names = F)
    message(paste("\n\n    *   files saved! old file removed   *\n\nfile name: 'table_for_metaboanalyst_",date,".csv'\n\n",
                  sep = ""))
    
  } else { 
    message(paste("\nerror:  the file 'table_for_metaboanalyst_",date,
                  ".csv' is already exist\n  *change the files name if you want to keep as backup", sep = ""))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
} else {
  write.csv(lc.for.met,
            paste(path.for.saving.files,"/table_for_metaboanalyst_",date,".csv", sep = ""),
            row.names = F)
  message(paste("\n\n\t*   files saved!   *\n\nfile name: 'table_for_metaboanalyst_",date,".csv'\n\n",
                sep = ""))
}
}
