merge_with_library <- function(table.from.lcms.software, path.for.saving.file,
                                  path.to.library.file = "P:/pomegranate RNAseq 2022/lcms files/pom_library_updated.xlsx")
  {

##############
  
## table.from.lcms.software
# raw data from lcms (excel file in csv format from lcms in Tel-Hai)
#   example: table.from.lcms.software = "P:/yonatan/RNAseq_yonatan_2021/lcms/data.from.lcms/All_Compounds_pom_T_060722.csv"


## path for the saved file directory
#   example: path.for.saving.file = "P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/transgenic/transgenic_files"

##############

# remove "/" if its in the end of "path.for.saving.file" row
  if (substr(path.for.saving.file, nchar(path.for.saving.file)-1+1, nchar(path.for.saving.file)) == "/") {
    path.for.saving.file = substr(path.for.saving.file,1,nchar(path.for.saving.file)-1)
}

date = gsub("/","",format(Sys.time(), "%x"))

lc.data = read.csv(table.from.lcms.software)
names(lc.data)[1] = "Tags"
lc.mzcl = lc.data[!lc.data$Tags == "",]
lc = data.frame(lc.mzcl$Name)
names(lc) = "Name"
lc$Name = gsub("std", "STD", lc$Name)
lc = unique(lc, by = "Name")

for (i in 1:length(lc[,1])) {
  lc[i,] = gsub("\r|\n","",lc[i,])
}



name2pathway = readxl::read_excel(path.to.library.file)
names(name2pathway)[1] = "Name"
name2pathway = unique(name2pathway, by = "Name")

for (i in 1:length(name2pathway$Name)) {
  name2pathway$Name[i] = gsub("\r|\n","",name2pathway$Name[i])
}


new.comp.and.pathway.list = merge.data.frame(lc, name2pathway, by = "Name", all.x = T)



if (file.exists(paste(path.for.saving.file,"/merge.comp.with.library.",date,".csv", sep = "")) == T) {
  message(paste("the file 'merge.comp.with.library.",date,
                ".csv' is already exist, do you want to continue?\n  it will remove the old file", sep = ""))
  
  if (regexpr(readline("   (y/n): "), 'y', ignore.case = TRUE) == 1) {
    write.csv(new.comp.and.pathway.list,
              paste(path.for.saving.file,"/merge.comp.with.library.",date,".csv", sep = ""),
              row.names = F, na = "")
    message(paste("\n\n*  file saved in your path. old file removed  *\n\nfile name: merge.comp.with.library.",date,".csv", sep = ""))
  } else { 
    message(paste("\nerror:  the file 'merge.comp.with.library.",date,
                  ".csv' is already exist\n\nchange the files names if you want to keep as backup\n", sep = ""))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
}  else {
  write.csv(new.comp.and.pathway.list,
            paste(path.for.saving.file,"/merge.comp.with.library.",date,".csv", sep = ""),
            row.names = F, na = "")
  message(paste("\n\n*  file saved in your path  *\n\nfile name: merge.comp.with.library.",date,".csv\n\n", sep = ""))
  }
}
