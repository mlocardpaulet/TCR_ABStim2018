################################################################################
# Correlation analysis of the Merged samples.
################################################################################



################################################################################
###### Import table
################################################################################

load("RData/08_StatResults_Phospho.RData")

################################################################################
###### Packages
################################################################################
library(reshape2)
library(Hmisc) # Pearson correlation
library(matrixStats) # rowMax()
library(ape) # as.phylo()
library(gplots)

################################################################################
# For this analysis I replace missing values with the mean values of 200 loops. See the Scripts in the folder "StatisticalAnalysis" for more information.
################################################################################

tab <- export[,grepl("MeanLoops", names(export))]
tabtech <- export[,grepl("MeanTech", names(export))]
row.names(tab) <- export$psiteID

# Keep only the sites that have a minimum of 5 measured values in a minimum of one data set (TiO2 or pY-IP).
k <- sapply(seq_len(nrow(tabtech)), function(x) {
  length(tabtech[x,][!is.na(tabtech[x,])])
})
tab <- tab[k >= 5,]

########################################################################
# I perform the correlation calculations.

# I keep only the regulated phosphorylation sites:
mat <- tab[row.names(tab) %in% as.character(export$psiteID[export$Regulation == "TRUE"]),]

mat2 <- t(scale(t(mat))) # rowwise scaling

# Combine the pY and the TiO2 data:
pY <- mat2[,grepl("pTyr", colnames(mat2))]
pY <- pY[sapply(seq_len(nrow(pY)), function(x) {length(pY[x,][!is.na(pY[x,])]) > 0}),]
TiO2 <- mat2[,grepl("TiO2", colnames(mat2))]
TiO2 <- TiO2[sapply(seq_len(nrow(TiO2)), function(x) {length(TiO2[x,][!is.na(TiO2[x,])]) > 0}),]
colnames(pY) <- gsub("_pTyr", "", colnames(pY), fixed = T)
colnames(TiO2) <- gsub("_TiO2", "", colnames(TiO2), fixed = T)
pY <- pY[,order(colnames(pY))]
TiO2 <- TiO2[,order(colnames(TiO2))]
pY <- cbind(pY[,1:6], matrix(NA, ncol = 6, nrow = nrow(pY)), pY[,7:18])
colnames(pY) <- colnames(TiO2)
mat2 <- rbind(pY, TiO2)

# Perform the pearson correlation between phosphorylation sites:
pearsCor <- rcorr(t(mat2), type="pearson")

matR <- pearsCor[[1]]
matP <- pearsCor[[3]]

save(list = c("matR", "matP"), file = "RData/09_Correlations.Rdata")

colfunclight <- colorRampPalette(c("darkblue", "blue", "deepskyblue", "white", "yellow", "red", "darkred"))

pdf("Figures/Correlation.pdf", useDingbats=FALSE, 14, 12)
heatmap.2(matR, trace = "n", col = colfunclight(100), cexRow = 0.25, cexCol = 0.25)
dev.off()


################################################################################
# Export for cytoscape
matR <- melt(matR)
matP <- melt(matP)
matRP <- cbind(matR, matP[,ncol(matP)])
names(matRP) <- c("Psite1", "Psite2", "R", "pvalue")
#matRP <- matRP[matRP$pvalue<=0.05,]
#matRP <- matRP[abs(matRP$R)>=0.9,]
matRP <- matRP[!is.na(matRP$Psite1),]

# For CytoscapeFiltering:
matRP$PvalUnder0.01 <- matRP$pvalue <= 0.01

matRP <- data.frame(matRP, "-log10pval" = -(log(matRP$pvalue,10)))

matRP <- matRP[matRP$PvalUnder0.01 != FALSE,]

# Remove duplicated connections:
matRP <- matRP[!is.na(matRP$R),]
for(el in as.character(unique(matRP$Psite1))) {
       for (el2 in as.character(unique(matRP$Psite2[matRP$Psite1==el]))) {
               matRP <- matRP[!(matRP$Psite1==el2 & matRP$Psite2==el),]
       }
}
matRP <- matRP[!duplicated(matRP),]

write.table(matRP[!is.na(matRP$R),], "SupTables/Cytoscape_KineticsCorrelation.txt", sep = "\t", row.names = F, quote = F)

