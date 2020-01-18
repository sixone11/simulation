#!/usr/bin/env Rscript
path <-"input" 
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){paste(path,x,sep='/')})  
data <- lapply(filePath, function(x){read.table(x, header=F,sep="\t")})  

for(i in 1:length(data)){
  data[[i]]<-data[[i]]$V5
}


data_process<-data[[1]]
for(i in 2:length(data)){
  data_process<-cbind(data_process,data[[i]])
}

####################################################################################
threhold<-0.5
data_new<-data_process
data_new[data_new>0.3]=1
data_new[data_new<=0.3]=0
data_new[is.na(data_new)]=0
correlation<-cor(data_new,method="kendall")
R2<-correlation^2
r_1<-1-R2
matrix_ken<-cmdscale(r_1, k=2)
for(i in 1:length(data)){
  rownames(matrix_ken)[i]<-strsplit(x=names(data)[i], split="_rate")[[1]][1]
}
MDS<-data.frame(matrix_ken)
MDS$software<-rownames(matrix_ken)
colnames(MDS)<-c("x","y","software")
write.table(MDS,file="output/kendallW_results",append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(cowplot))
suppressMessages(library(ggthemes))
suppressMessages(library(grid))

p<-ggplot(MDS,aes(x=x,y=y,color=software,color=software,shape=software,stroke=1.5),alpha = 0.5)+geom_point(size=4,stroke=2)+labs(x="1st coordinate",y="2nd coordinate")+theme(axis.text.x=element_text(size = 12,face = "bold"),axis.title=element_text(face="bold"),legend.title = element_text(face="bold"),axis.text.y=element_text(size=12,face = "bold"),axis.line=element_blank())+theme_bw()+theme(axis.title =element_text(size = 24),axis.text =element_text(size = 20, color = 'black'))+labs(color="methods")+scale_colour_manual(name = "methods",labels = c("BiSeq","DMRfinder", "DSS", "HMM-DM","HMM-Fisher","methylKit","methylSig","metilene"),values = c(rgb(165,84,44,max=255),rgb(250,139,101,max=255),rgb(143,161,201,max=255),rgb(231,137,193,max=255),rgb(164,216,95,max=255),rgb(253,216,68,max=255),rgb(228,195,151,max=255),rgb(179,179,179,max=255)))+scale_shape_manual(name = "methods",labels = c("BiSeq","DMRfinder", "DSS", "HMM-DM","HMM-Fisher","methylKit","methylSig","metilene"),values=c(9,6,8,11,3,4,23,24))

pdf("MDS_different_threhold.pdf",width=18,height=12)
plot(p)
dev.off()