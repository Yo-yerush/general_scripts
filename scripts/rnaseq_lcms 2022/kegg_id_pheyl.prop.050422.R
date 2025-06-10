setwd("C:/Users/yonye/OneDrive/שולחן העבודה/מעבדה/RNAseq_yonatan_2021/lcms/lcms test 310322/")
lc.data = read.csv("Compounds_KeggID_pom_T_310322.csv")


lc.k = lc.data[lc.data[,5] == "Biosynthesis of phenylpropanoids",] #

data.by.row = data.frame(rownames(lc.k),lc.k[,c(9,10)])

write.csv(file.keggID, "keggID.transgenic.pre.310322.csv", row.names = F)


