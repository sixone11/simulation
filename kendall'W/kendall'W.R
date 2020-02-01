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
threshold<-0.5
data_new<-data_process
data_new[data_new>threshold]=1
data_new[data_new<=threshold]=0
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
rownames(correlation)<-rownames(matrix_ken)
colnames(correlation)<-rownames(matrix_ken)
write.table(correlation,file="output/kendallW_results",quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
write.table(MDS,file="output/MDS_results",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
