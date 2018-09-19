# load(file = "RData/12_PhosphoTableWithProteinWaring.RData")
# PSP <- read.table(file = "SupTables/PhosphoSitePlusMapping.txt", header = T, sep = "\t", stringsAsFactors = F)
# tab <- export
# tab$Kinase <- PSP$Kinase[match(tab$psiteID, PSP$ID)]

coloursHM <- colorRampPalette(c("green", "black", "red"))
require(reshape2)
require(gplots)

PlotHM <- function(psites, tab, save = T) { # psites is the vector of phosphosite ID to plot, tab is the table with the appropriate quantification values
  cols <- which(grepl("MeanLoops", names(tab)))
  gtab <- tab[tab$psiteID %in% psites,cols]
  row.names(gtab) <- tab$GeneID[tab$psiteID %in% psites]
  colnames(gtab) <- gsub("MeanLoops_", "", colnames(gtab), fixed = T)
  pY <- gtab[,grepl("pTyr", colnames(gtab))]
  TiO2 <- gtab[,grepl("TiO2", colnames(gtab))]
  kpY <- sapply(seq_len(nrow(pY)), function(x) {
    length(pY[x,][!is.na(pY[x,])]) > 0
  })  
  kTiO2 <- sapply(seq_len(nrow(TiO2)), function(x) {
    length(TiO2[x,][!is.na(TiO2[x,])]) > 0
  })  
  pY <- pY[kpY,]
  TiO2 <- TiO2[kTiO2,]
  ltab <- list(pY, TiO2)
  for (el in ltab) {
    if (nrow(el) > 0) {
      na <- row.names(el)
      el <- apply(el, 2, as.numeric)
      row.names(el) <- na
      regulation <- tab$Regulation[match(na, tab$GeneID)]
      if (sum(grepl("TiO2", colnames(el))) > 0) {
        colnames(el) <- gsub("R[12345]_TiO2_", "", colnames(el))
        colnames(el) <- gsub("S", "", colnames(el))
        colnames(el) <- gsub(".", "", colnames(el), fixed = T)
        colnames(el)[colnames(el)=="N"] <- 0
        el <- el[,order(as.numeric(colnames(el)))]
        heatmap.2(el, col = coloursHM(100), na.color = "grey30", colsep = c(4, 8, 12, 16, 20), trace = "n", margins = c(8,8), Colv = F)
        if (save == T) {
          na <- paste0("Figures/Heatmaps/", Sys.time(), ".jpg")
          na <- gsub(":", "", na, fixed = T)
          pdf(na, useDingbats=FALSE, height = 11.69, width = 8.27)
          heatmap.2(el, col = coloursHM(100), na.color = "grey30", colsep = c(4, 8, 12, 16, 20), trace = "n", margins = c(8,8))
          dev.off()
        }
      } else {
        colnames(el) <- gsub("R[12345]_pTyr_", "", colnames(el))
        colnames(el) <- gsub("S", "", colnames(el))
        colnames(el) <- gsub(".", "", colnames(el), fixed = T)
        colnames(el)[colnames(el)=="N"] <- 0
        el <- el[,order(as.numeric(colnames(el)))]
        if (save == T) {
          na <- paste0("Figures/Heatmaps/", Sys.time(), ".jpg")
          na <- gsub(":", "", na, fixed = T)
          pdf(na, useDingbats=FALSE, height = 11.69, width = 8.27)
          heatmap.2(el, col = coloursHM(100), na.color = "grey30", colsep = c(3, 6, 9, 12), trace = "n", margins = c(8,8))
          dev.off()
        }
      }  
    }
  }
}


################################################################################

coloursHMFC <- colorRampPalette(c("darkblue", "white","darkred"))

PlotHMmeans <- function(psites, tab, save = T) { # psites is the vector of phosphosite ID to plot, tab is the table with the appropriate quantification values
  cols <- which(grepl("T0.FC", names(tab), fixed = T))
  gtab <- tab[tab$psiteID %in% psites,cols]
  row.names(gtab) <- tab$GeneID[tab$psiteID %in% psites]
  colnames(gtab) <- gsub("T", "", colnames(gtab), fixed = T)
  colnames(gtab) <- sapply(colnames(gtab), function(x) {
    strsplit(x, ".", fixed = T)[[1]][1]
  })
  gtab <-cbind(gtab, "0" = rep(0, nrow(gtab)))
  gtab <- gtab[,order(as.numeric(colnames(gtab)))]
  regulation <- tab$Regulation[match(row.names(gtab), tab$GeneID)]
  regulation <- ifelse(regulation == "TRUE", 2, 1)
  na <- row.names(gtab)
  gtab <- apply(gtab, 2, as.numeric)
  row.names(gtab) <- na
  heatmap.2(gtab, col = coloursHMFC(100), na.color = "grey30",trace = "n", margins = c(8,8), Colv = F, colRow = regulation)
  if (save == T) {
    na <- paste0("Figures/Heatmaps/", Sys.time(), ".pdf")
    na <- gsub(":", "", na, fixed = T)
    pdf(na, useDingbats=FALSE, height = 11.69, width = 8.27)
    heatmap.2(gtab, col = coloursHMFC(100), na.color = "grey30",trace = "n", margins = c(8,8), Colv = F, colRow = regulation)
    dev.off()
  }
}

