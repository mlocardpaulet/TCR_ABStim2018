################################################################################
#      20180926 - export sequence windows for motif enrichment analysis        #
################################################################################

load("T:/TCR/LabelFreeCD3CD4/TCR_ABStim2018/RData/05_ParsedTablesNormPSP.RData")

tab <- lapply(lf, function(x) {
  x[,c("phosphoSites","Sequence.window")]
})
source(file = "RFunctions/RBindList.R")
tab <- RBindList(tab)
tab <- tab[!duplicated(tab),]
tab$Sequence.window7 <- as.character(sapply(tab$Sequence.window, substr, 9, 23))

write.table(tab, file = "SupTables/SequenceWindows.txt", sep = "\t", row.names = F, quote = F)
