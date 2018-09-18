# load(file = "RData/08_StatResults_Phospho.RData")
# tab <- export
# tab <- read.table(file = "SupTables/PhosphoData.txt", sep = "\t", header = T)
coloursLines <- c("blue4", "deepskyblue4", "lightseagreen", "mediumspringgreen")
require(ggplot2)
require(reshape2)

PlotSite <- function(psites, tab, save = T) { # psites is the vector of phosphosite ID to plot, tab is the table with the appropriate quantification values
  cols <- which(grepl("MeanLoops", names(tab)))
  for (el in psites) {
    gtab <- tab[tab$psiteID==el,cols]
    geneID <- unique(tab$GeneID[tab$psiteID == el])
    names(gtab) <- gsub("MeanLoops_", "", names(gtab), fixed = T)
    gtab <- melt(gtab)
    if (length(gtab$value[!is.na(gtab$value)]) > 0) {
      gtab$variable <- as.character(gtab$variable)
      gtab$TimePoint <- sapply(gtab$variable, function(x) {
        strsplit(x, "_", fixed = T)[[1]][3]
      })
      gtab$TimePoint[grepl("NS.", gtab$TimePoint)] <- 0 
      gtab$TimePoint[grepl("S30.", gtab$TimePoint, fixed = T)] <- 30 
      gtab$TimePoint[grepl("S15.", gtab$TimePoint, fixed = T)] <- 15 
      gtab$TimePoint[grepl("S120.", gtab$TimePoint, fixed = T)] <- 120 
      gtab$TimePoint[grepl("S300.", gtab$TimePoint, fixed = T)] <- 300 
      gtab$TimePoint[grepl("S600", gtab$TimePoint, fixed = T)] <- 600 
      gtab$Replicate <- sapply(gtab$variable, function(x) {
        strsplit(x, "_", fixed = T)[[1]][1]
      })
      gtab$TimePoint <- as.numeric(as.character(gtab$TimePoint))
      gtab <- gtab[!is.na(gtab$value),]
      
      g <- ggplot(gtab, aes(x = TimePoint, y = value)) + geom_line(aes(x = TimePoint, y = value, group = Replicate, col = Replicate), size = 1.2) + geom_vline(xintercept = 0, size = 1.2) + geom_point(col = "black", shape = "+", size = 5, alpha = 0.8) + theme_minimal() + ggtitle(paste0(el, " ", geneID)) + scale_color_manual(values = coloursLines[seq_along(unique(gtab$Replicate))]) + ylab("log2-transformed normalised MS intensities") # + geom_boxplot(aes(x = TimePoint, y = value), width = 0.5)
    } else {
      print(paste0("No quantification values for ", el, "."))
    }
    print(g)
    if (save == T) {
      na <- paste0("Figures/Kinetics/", el, ".jpg")
      na <- gsub("|", "", na, fixed = T)
      jpeg(na, 5.845, 4.135, units = "in", res = 300) # width and height in inches.
      print(g)
      dev.off()
    }
  }
}

################################################################################

PlotSiteYLim <- function(psites, tab, ylim = c(10,40)) { # psites is the vector of phosphosite ID to plot, tab is the table with the appropriate quantification values
  # Here I can input the ylim to have it fixed.
  # The plot is printed but not saved
  cols <- which(grepl("MeanLoops", names(tab)))
  for (el in psites) {
    gtab <- tab[tab$psiteID==el,cols]
    geneID <- unique(tab$GeneID[tab$psiteID == el])
    names(gtab) <- gsub("MeanLoops_", "", names(gtab), fixed = T)
    gtab <- melt(gtab)
    if (length(gtab$value[!is.na(gtab$value)]) > 0) {
      gtab$variable <- as.character(gtab$variable)
      gtab$TimePoint <- sapply(gtab$variable, function(x) {
        strsplit(x, "_", fixed = T)[[1]][3]
      })
      gtab$TimePoint[grepl("NS.", gtab$TimePoint)] <- 0 
      gtab$TimePoint[grepl("S30.", gtab$TimePoint, fixed = T)] <- 30 
      gtab$TimePoint[grepl("S15.", gtab$TimePoint, fixed = T)] <- 15 
      gtab$TimePoint[grepl("S120.", gtab$TimePoint, fixed = T)] <- 120 
      gtab$TimePoint[grepl("S300.", gtab$TimePoint, fixed = T)] <- 300 
      gtab$TimePoint[grepl("S600", gtab$TimePoint, fixed = T)] <- 600 
      gtab$Replicate <- sapply(gtab$variable, function(x) {
        strsplit(x, "_", fixed = T)[[1]][1]
      })
      gtab$TimePoint <- as.numeric(as.character(gtab$TimePoint))
      gtab <- gtab[!is.na(gtab$value),]
      
      g <- ggplot(gtab, aes(x = TimePoint, y = value)) + geom_line(aes(x = TimePoint, y = value, group = Replicate, col = Replicate), size = 1.2) + geom_vline(xintercept = 0, size = 1.2) + geom_point(col = "black", shape = "+", size = 5, alpha = 0.8) + theme_minimal() + ggtitle(paste0(el, " ", geneID)) + scale_color_manual(values = coloursLines[seq_along(unique(gtab$Replicate))]) + ylab("log2-transformed normalised MS intensities") + ylim(ylim) # + geom_boxplot(aes(x = TimePoint, y = value), width = 0.5)
    } else {
      print(paste0("No quantification values for ", el, "."))
    }
    print(g)
  }
}

################################################################################

# psites <- as.character(tab$psiteID[tab$Regulation == "TRUE"][order(tab$pAnova[tab$Regulation == "TRUE"], decreasing = T)][1:30])
# psites <- as.character(tab$psiteID[tab$ClusterMerged == "9"])
# psites <- psites[!is.na(psites)]
# psites <- c("Q3TTA7_Y889")
# psites <- as.character(tab$psiteID[grepl("Q60949", tab$psiteID)])
# keep <- which(tab$Regulation == "TRUE" & tab$ProportioPassStat < 1)
# psites <- as.character(tab$psiteID[keep])
# PlotSite(psites, tab)
# PlotSite(tab$psiteID[tab$GeneID %in% c("Cd5_Y452")], tab)
# PlotSite(c("Q3V3E1_Y9"), tab)
# PlotSite(tab$psiteID[grepl("Ubash", tab$GeneID)], tab)
