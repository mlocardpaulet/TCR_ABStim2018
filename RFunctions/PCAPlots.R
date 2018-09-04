# Visualisation of PCA results

plotSDev <- function(pca) {
  tab <- unclass(summary(pca))$importance
  CumProp <- tab[3,]
  xlab <- paste0(rep("Dim ", length(CumProp)), seq_along(CumProp))
  barplot(CumProp, main = "Cumulative proportion of variance per dimention", names.arg = xlab)
}
  
plotConditions <- function(pca, dim1 = 1, dim2 = 2) {
  gtab <- as.data.frame(unclass(pca$rotation), stringsAsFactor = F)
  gtab$Condition <- row.names(gtab)
  print(ggplot(gtab, aes(x = gtab[,dim1], y = gtab[,dim2], label = Condition)) + theme_bw()+ geom_point(size = 3) + ylab(paste0("Dim. ", dim1, "(", summary(pca)$importance[2,dim1]*100, "%)")) + xlab(paste0("Dim. ", dim2, "(", summary(pca)$importance[2,dim2]*100, "%)"))  + geom_text_repel())
}

plotGenes <- function(pca, dim1 = 1, dim2 = 2) {
  gtab <- as.data.frame(unclass(pca$x), stringsAsFactor = F)
  gtab$Gene <- row.names(gtab)
  print(ggplot(gtab, aes(x = gtab[,dim1], y = gtab[,dim2])) + theme_bw()+ geom_point(size = 3, alpha = 0.3) + ylab(paste0("Dim. ", dim1, "(", summary(pca)$importance[2,dim1]*100, "%)")) + xlab(paste0("Dim. ", dim2, "(", summary(pca)$importance[2,dim2]*100, "%)")))
  smoothScatter(gtab[,dim1], gtab[,dim2], nrpoints = Inf, cex = 2)
}

# 
# # Manual:
# gtab <- as.data.frame(unclass(pca$x), stringsAsFactor = F)
# gtab$Gene <- row.names(gtab)
# gtab$Gene[!grepl("Y", gtab$Gene)] <- NA
# ggplot(gtab, aes(x = gtab[,dim1], y = gtab[,dim2], label = Gene)) + theme_bw()+ geom_point(size = 3, alpha = 0.3) + ylab(paste0("Dim. ", dim1, "(", summary(pca)$importance[2,dim1]*100, "%)")) + xlab(paste0("Dim. ", dim2, "(", summary(pca)$importance[2,dim2]*100, "%)")) + geom_text_repel()
# smoothScatter(gtab[,dim1], gtab[,dim2], nrpoints = Inf, cex = 2)