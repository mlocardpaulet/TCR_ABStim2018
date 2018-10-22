#plot kinetic of cluster by keywords
library(xlsx)
library(lattice)
library(ggplot2)
library("RColorBrewer")
library(reshape)
library(devtools)
library(gridExtra)
library(stringr)
library(svglite)
library(Rtsne)
library(ggrepel)
library(plotly)
library(ClusterX)


#mettre "kineticPsite1.csv" pour tous les sites

setwd("S:/Romain/IN PROGRESS/SYBILLA/phosphoproteome/analyse Phosphoproteome anti-CD3-CD4/2017/clusters")
kinetic<-read.csv("kineticPsiteNew2018.csv",sep=";", dec=",",stringsAsFactors=FALSE)

setwd("S:/Romain/IN PROGRESS/SYBILLA/phosphoproteome/analyse Phosphoproteome anti-CD3-CD4/2017/R figurgene")
phositeprot<-read.csv("liste phosphosite-proteinsNew.csv",sep=";", dec=",",stringsAsFactors=FALSE)
phositeprot$Gene<-toupper(phositeprot$Gene)

phositeprot<-subset(phositeprot,phositeprot$Regulated=="TRUE")
kinetic<-subset(kinetic,kinetic$Regulation =="TRUE")

t0<-grep("NS",colnames(kinetic))
t15<-grep("S15",colnames(kinetic))
t120<-grep("S120",colnames(kinetic))
t300<-grep("S300",colnames(kinetic))
t600<-grep("S600",colnames(kinetic))
t30<-grep("S30",colnames(kinetic))
t30<-setdiff(t30,t300)

kineticMean<-data.frame(matrix(NA,length(kinetic[,1]),9))
colnames(kineticMean)<-c("psiteID","NS","15","30","120","300","600","cluster","gene")

kineticMean$psiteID<-kinetic$psiteID
kineticMean$NS<-rowMeans(kinetic[,t0],na.rm=TRUE)
kineticMean$`15`<-rowMeans(kinetic[,t15],na.rm=TRUE)
kineticMean$`30`<-rowMeans(kinetic[,t30],na.rm=TRUE)
kineticMean$`120`<-rowMeans(kinetic[,t120],na.rm=TRUE)
kineticMean$`300`<-rowMeans(kinetic[,t300],na.rm=TRUE)
kineticMean$`600`<-rowMeans(kinetic[,t600],na.rm=TRUE)

X<-1
for (i in 1:length(kinetic[,1])){
  
  n<-which(kineticMean[X,1]==phositeprot$psiteID)
  if (length(n)==1)  {kineticMean$cluster[X]<-phositeprot$Cluster[n] 
  kineticMean$gene[X]<-phositeprot$Gene[n] 
  }
  else kineticMean$cluster[X]<--1
  
  X<-X+1}

kineticMeanN<-data.frame(matrix(NA,length(kinetic[,1]),9))
colnames(kineticMeanN)<-c("psiteID","NS","15","30","120","300","600","cluster","gene")
kineticMeanN$psiteID<-kineticMean$psiteID
kineticMeanN$cluster<-kineticMean$cluster
kineticMeanN$gene<-kineticMean$gene


kineticMean[,2:7]<-2^kineticMean[,2:7]

X<-1
for (i in 1:length(kinetic[,1])){
  
  #mean or max or z-score
  moyp<-mean(as.numeric(kineticMean[X,2:7]))
  maxp<-max(kineticMean[X,2:7])
  st<-sd(as.numeric(kineticMean[X,2:7]))
  #here
  #kineticMeanN[X,2:7]<-kineticMean[X,2:7]/maxp
  #kineticMeanN[X,2:7]<-kineticMean[X,2:7]/moyp 
  kineticMeanN[X,2:7]<-(kineticMean[X,2:7]-moyp)/st
  
  X<-X+1}

#modify double values 
kineticMeanN[256,2]<-0.37
kineticMeanN[256,2]<-0.93


a<-Rtsne(kineticMeanN[,2:7],initial_dims = 6)
b<-as.data.frame(a$Y)
b<-cbind(kineticMean$psiteID,b,kineticMean$cluster,kineticMean$gene)
b[b==-1]<-10
b$`kineticMean$psiteID`<-gsub(";", "|", b$`kineticMean$psiteID`)


#setwd("S:/Romain/IN PROGRESS/SYBILLA/phosphoproteome/analyse Phosphoproteome anti-CD3-CD4/2017/t-SNE")
#write.xlsx(b,"table_all_phospho_loca_tSNENew2018 ClusterX.xlsx")
#si recuperation all coordo
setwd("S:/Romain/IN PROGRESS/SYBILLA/phosphoproteome/analyse Phosphoproteome anti-CD3-CD4/2017/t-SNE")
b<-read.csv("table_all_phospho_loca_tSNENew2018.csv",sep=";", dec=",",stringsAsFactors=FALSE)
b<-b[-1]
colnames(b)<-c("kineticMean$psiteID","V1","V2","kineticMean$cluster","kineticMean$gene")


# ClusterX
d<-ClusterX(b[,2:3])

l<-data.frame(matrix(NA,length(b[,1]),5))
colnames(l)<-colnames(b)
l<-b
l$`kineticMean$cluster`<-d$cluster


b$`kineticMean$cluster`<-l$`kineticMean$cluster`


cc<-ggplot(b,aes(b$V1,b$V2,color=factor(b$`kineticMean$cluster`),label=b$`kineticMean$psiteID`),colour=factor(b$`kineticMean$cluster`))+
  geom_point(size=5) +
  
  
  ggtitle("t-SNE") +
  theme(legend.position="none")

cc

