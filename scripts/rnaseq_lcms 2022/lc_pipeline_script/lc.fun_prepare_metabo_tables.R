prepare_metabo_tables <- function(table.from.lcms.software, library.merged.file,
                                  lcms.samples.file, path.for.saving.files, 
                                  first_row_phenotype, sample_order_name = F,
                                  Tags = T) {

  ##############
  
  ### table.from.lcms.software
  # full path to raw data file from lcms (csv file from lcms software in Tel-Hai)
  #   example: table.from.lcms.software = "P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/All_Compounds_pom_P_060722.csv"
  
  
  ### library.merged.file
  # full path to the merged library file you created last step with the "lc_merge_with_library" function
  
  ### lcms.samples.file
  # full path to samples file (csv) in order of lcms running, with BLK and QC
  #   example: lcms.samples.file = "P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/LCMS_rnaseq_samples/LCMS_rnaseq.samples_190122.csv"
  
  
  ### path.for.saving.files 
  # path for the saved files directory
  #   example: path.for.saving.files = "P:/yonatan/RNAseq_yonatan_2021/lcms/csv_files_for_metabo/comp.with.pathway.P/" 
  
  
  ### first_row_phenotype
  # group description for metaboanalyst analysis
  #   example for 30 samples and 3 groups: 
  # first_row_phenotype = c("OE_plants", rep("wild.type",10), rep("mt.line1",10), rep("mt.line2",10))
  
  ### sample_order_name 
  # to get sample names easily, run script with: 
  #   all  non-default arguments except <first_row_phenotype> and with <sample_order_name = T>
  #   *no file will saved
  
  
  ### Tags
  # logical, if used "Tags" column for filtering at the lcms software in Tel-Hai
  
  ##############
  
  date = gsub("/","",format(Sys.time(), "%x"))
  
# remove "/" if its in the end of "path.for.saving.files" row
if (substr(path.for.saving.files, nchar(path.for.saving.files)-1+1, nchar(path.for.saving.files)) == "/") {
  path.for.saving.files = substr(path.for.saving.files,1,nchar(path.for.saving.files)-1)
}


df.1 = read.csv(library.merged.file)[,1:2]
names(df.1)[1] = "Name"
df.1 = df.1[!df.1$comp.pathway == "",]
for (i in 1:length(df.1$Name)) {
  df.1$Name[i] = gsub("\r|\n","",df.1$Name[i])
}

samples = read.csv(lcms.samples.file)

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

# merge to create table for.metabo and for names.with.pathway
comp.pathway.name.table.with.metabo = merge.data.frame(lc.for.met,df.1, by = "Name", all.y = T)

# add puni.marker to the 1st row
comp.pathway.name.table.with.metabo[(nrow(comp.pathway.name.table.with.metabo)+1),1:(length(comp.pathway.name.table.with.metabo)-1)] = first_row_phenotype
comp.pathway.name.table.with.metabo = comp.pathway.name.table.with.metabo[c(nrow(comp.pathway.name.table.with.metabo),
                                                                            1:(nrow(comp.pathway.name.table.with.metabo)-1)),]

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

# make non-duplicates "comp.pathway" by adding 'X[d]'
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

if (file.exists(paste(path.for.saving.files,"/files_for_metaboanalyst_",date, sep = "")) == T) {
  message(paste("the folder 'files_for_metaboanalyst_",date,
                             "' is already exist, do you want to continue?\n  *it will remove old files with the same name in that folder", sep = ""))
  
  if (regexpr(readline("\t(y/n): "), 'y', ignore.case = TRUE) == 1) {
    dir.create(paste(path.for.saving.files,"/files_for_metaboanalyst_",date, sep = ""))
  } else { 
    message(paste("\nerror:  the folder 'files_for_metaboanalyst_",date,
                  "' is already exist\n  *change the files/folder names if you want to keep as backup", sep = ""))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
}  else { dir.create(paste(path.for.saving.files,"/files_for_metaboanalyst_",date, sep = ""))}

  write.csv(file1,
            paste(path.for.saving.files,"/files_for_metaboanalyst_",date,"/pathway.name.for.metabo.",date,".csv", sep = ""),
            row.names = F)
  write.csv(comp.pathway.name.table.with.metabo[,1:(length(comp.pathway.name.table.with.metabo)-1)],
            paste(path.for.saving.files,"/files_for_metaboanalyst_",date,"/comp.name.for.metabo.",date,".csv", sep = ""),
            row.names = F)
  write.csv(comp.pathway.name.table.with.metabo[-1,c(1,length(comp.pathway.name.table.with.metabo))],
            paste(path.for.saving.files,"/files_for_metaboanalyst_",date,"/comp.and.pathway.table.",date,".csv", sep = ""),
            row.names = F)
  
  message(paste("\n\n   *  files saved in 'files_for_metaboanalyst_",date,"' folder  *\n\n",
                      sep = ""))
}
