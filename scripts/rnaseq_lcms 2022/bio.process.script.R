yo1 = read.csv("../../DESeq2/transgenic/statistics/high_exp_bHLH94_VS_EV/high_exp_bHLH94_VS_EV.FC.up.DE.csv")
yo3 = yo1$biological.process
write.csv(yo3,"../../DESeq2/transgenic/statistics/high_exp_bHLH94_VS_EV/bio.process.vector.pre.high_exp_bHLH94_VS_EV.FC.up.DE.csv",row.names = F)
yo1.2 = read.csv("../../DESeq2/transgenic/statistics/high_exp_bHLH94_VS_EV/high_exp_bHLH94_VS_EV.FC.down.DE.csv")
yo3.2 = yo1.2$biological.process
write.csv(yo3.2,"../../DESeq2/transgenic/statistics/high_exp_bHLH94_VS_EV/bio.process.vector.pre.high_exp_bHLH94_VS_EV.FC.down.DE.csv",row.names = F)




yo4 = read.csv("../../DESeq2/transgenic/statistics/bio.process.high_exp_bHLH94_VS_EV.FC.up.DE.csv")
yo5 = read.csv("../../DESeq2/transgenic/statistics/bio.process.high_exp_bHLH94_VS_EV.FC.down.DE.csv")

yo4.2 = table(yo4)
yo5.2 = table(yo5)

length(yo5.2)=length(yo4.2)

yo6 = data.frame(yo4.2,yo5.2)
yo6 = yo6[,c(2,1,3)]
names(yo6) = c("up.file","biological.process","down.file")

write.csv(yo6,"../../DESeq2/transgenic/statistics/bio.process.teble.high_exp_bHLH94_VS_EV.csv",row.names = F)

