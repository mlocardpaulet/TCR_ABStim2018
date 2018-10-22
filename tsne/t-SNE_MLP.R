#plot kinetic of cluster by keywords
library(lattice)
library(ggplot2)
library("RColorBrewer")
library(reshape2)
library(devtools)
library(gridExtra)
library(stringr)
library(svglite)
library(Rtsne)
library(ggrepel)
library(plotly)

# if(!require(devtools)){
#   install.packages("devtools") # If not already installed
# }
# devtools::install_github("JinmiaoChenLab/ClusterX")

library(ClusterX)


load("RData/12_PhosphoTableWithProteinWaring.RData")
kinetic <- export[,grepl("MeanLoops", names(export))]
names(kinetic) <- gsub("MeanLoops_", "", names(kinetic), fixed = T)
row.names(kinetic) <- export$psiteID

# Merge the TiO2 and pY tables:
tabpY <- kinetic[,grepl("pTyr", names(kinetic))]
tabTiO2 <- kinetic[,grepl("TiO2", names(kinetic))]
names(tabpY) <-  gsub("_pTyr", "", names(tabpY))
names(tabTiO2) <-  gsub("_TiO2", "", names(tabTiO2))
temp <- matrix(ncol = 6, nrow = nrow(tabpY))
colnames(temp) <- gsub("R1", "R3", names(tabpY)[1:6])
tabpY <- cbind(tabpY, temp)
tabpY <- tabpY[,order(names(tabpY))]
rm(temp)
# Remove only missing values:
k <- sapply(seq_len(nrow(tabpY)), function(x) {
  length(tabpY[x,][is.na(tabpY[x,])]) < ncol(tabpY)
})
tabpY <- tabpY[k,]
k <- sapply(seq_len(nrow(tabTiO2)), function(x) {
  length(tabTiO2[x,][is.na(tabTiO2[x,])]) < ncol(tabTiO2)
})
tabTiO2 <- tabTiO2[k,]

kinetic <- rbind(tabpY, tabTiO2)

# Keep only the regulated sites:
kinetic <- kinetic[row.names(kinetic) %in% export$psiteID[export$Regulation == "TRUE"],]

t0<-grep("NS.",colnames(kinetic), fixed = T)
t15<-grep("S15.",colnames(kinetic), fixed = T)
t120<-grep("S120.",colnames(kinetic), fixed = T)
t300<-grep("S300.",colnames(kinetic), fixed = T)
t600<-grep("S600.",colnames(kinetic), fixed = T)
t30<-grep("S30.",colnames(kinetic), fixed = T)

kineticMean <- data.frame(matrix(NA,length(kinetic[,1]),7))
# colnames(kineticMean) <- c("psiteID","NS","15","30","120","300","600","cluster","gene")
names(kineticMean) <- c("psiteID","NS","15","30","120","300","600")

kineticMean$psiteID <- row.names(kinetic)
kineticMean$NS <- rowMeans(kinetic[,t0],na.rm=TRUE)
kineticMean$`15` <- rowMeans(kinetic[,t15],na.rm=TRUE)
kineticMean$`30` <- rowMeans(kinetic[,t30],na.rm=TRUE)
kineticMean$`120` <- rowMeans(kinetic[,t120],na.rm=TRUE)
kineticMean$`300` <- rowMeans(kinetic[,t300],na.rm=TRUE)
kineticMean$`600` <- rowMeans(kinetic[,t600],na.rm=TRUE)

# X <- 1
# for (i in 1:length(kinetic[,1])){
#   
#   n <- which(kineticMean[X,1] == phositeprot$psiteID)
#   if (length(n)==1)  {kineticMean$cluster[X]<-phositeprot$Cluster[n] 
#   kineticMean$gene[X]<-phositeprot$Gene[n] 
#   }
#   else kineticMean$cluster[X]<--1
#   
#   X<-X+1
# }

# kineticMeanN<-data.frame(matrix(NA,length(kinetic[,1]),9))
# colnames(kineticMeanN)<-c("psiteID","NS","15","30","120","300","600","cluster","gene")
# kineticMeanN$psiteID<-kineticMean$psiteID
# kineticMeanN$cluster<-kineticMean$cluster
# kineticMeanN$gene<-kineticMean$gene


kineticMean[,2:7] <- 2^kineticMean[,2:7]
kineticMeanScaled <- kineticMean
kineticMeanScaled[,2:7] <- t(scale(t(kineticMeanScaled[,2:7])))

kineticMeanN <- kineticMean
for (X in 1:length(kinetic[,1])) {
  #mean or max or z-score
  moyp <- mean(as.numeric(kineticMean[X,2:7]))
  maxp <- max(kineticMean[X,2:7])
  st <- sd(as.numeric(kineticMean[X,2:7]))
  #here
  #kineticMeanN[X,2:7]<-kineticMean[X,2:7]/maxp
  #kineticMeanN[X,2:7]<-kineticMean[X,2:7]/moyp 
  kineticMeanN[X,2:7] <- (kineticMean[X,2:7]-moyp)/st
}

# #modify double values 
# kineticMeanN[256,2]<-0.37
# kineticMeanN[256,2]<-0.93
# I artificially modify a bit the values of the site that is represented by 2 identical kinetics
kineticMeanN[kineticMeanN$psiteID == "Q8K4J6_S606",2:7] <- kineticMeanN[kineticMeanN$psiteID == "Q8K4J6_S606",2:7] + rnorm(n = 6, mean = 0.0001, sd = 0.001)

set.seed(546)
a <- Rtsne(kineticMeanN[,2:7],initial_dims = 6)
b <- as.data.frame(a$Y)
# b <- cbind(kineticMean$psiteID,b,kineticMean$cluster,kineticMean$gene)
b <- cbind(kineticMean$psiteID,b)
# b[b==-1] <- 10
# b$`kineticMean$psiteID` <- gsub(";", "|", b$`kineticMean$psiteID`)


#setwd("S:/Romain/IN PROGRESS/SYBILLA/phosphoproteome/analyse Phosphoproteome anti-CD3-CD4/2017/t-SNE")
#write.xlsx(b,"table_all_phospho_loca_tSNENew2018 ClusterX.xlsx")
#si recuperation all coordo
# setwd("S:/Romain/IN PROGRESS/SYBILLA/phosphoproteome/analyse Phosphoproteome anti-CD3-CD4/2017/t-SNE")
# b<-read.csv("table_all_phospho_loca_tSNENew2018.csv",sep=";", dec=",",stringsAsFactors=FALSE)
# b<-b[-1]
# colnames(b)<-c("kineticMean$psiteID","V1","V2","kineticMean$cluster","kineticMean$gene")


# ClusterX
d <- ClusterX(b[,2:3])

clusterPlot(d)
peakPlot(d)
densityPlot(d)


b$`kineticMean$cluster`<- d$cluster

cc <- ggplot(b, aes(b$V1, b$V2, color=factor(b$`kineticMean$cluster`), label=b$`kineticMean$cluster`), colour=factor(b$`kineticMean$cluster`)) +
  geom_point(size=5, alpha = 0.5) +
  ggtitle("t-SNE") +
  theme(legend.position="none") +
  theme_minimal() +
  geom_text()

cc

# reordering of the clusters to follow the kinetic:

reordering <- data.frame("old" = 1:15, "new" = c(1,6,13,8,11,5,12,2,7,14,15,10,4,9,3))
b$`kineticMean$cluster` <- reordering$new[match(b$`kineticMean$cluster`, reordering$old)]

cc <- ggplot(b, aes(b$V1, b$V2, color=factor(b$`kineticMean$cluster`), label=b$`kineticMean$psiteID`), colour=factor(b$`kineticMean$cluster`)) +
  geom_point(size=5, alpha = 0.5) +
  ggtitle("t-SNE") +
  theme(legend.position="none") +
  theme_minimal()

cc
