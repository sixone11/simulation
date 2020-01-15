args=commandArgs(T)
suppressMessages(library(ggplot2))
suppressMessages(library(pROC))
suppressMessages(library(plotROC))

calculate_auc <- function(name){
    file<-paste0(name,"_rate.bed")
    data<-read.table(file,header=F,sep="\t")
    AUC<-auc(data$V4, data$V5)
    p<-roc(data$V4,data$V5)
    plot_name<-paste0(name,".pdf")
    pdf("DMRfinder.pdf",width=6,height=6)
    plot(p)
    dev.off()
    return(data.frame(name=name,auc=AUC))
}
####################################################################
#debug


name=c("metilene","DMRfinder")
#"HMM-DM","HMM-Fisher","methylSig","methylKit","BiSeq","DSS")
results<-data.frame()
for(i in name){
  results<-rbind(results,calculate_auc(i))
}

write.table(results,file=args[1],append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


