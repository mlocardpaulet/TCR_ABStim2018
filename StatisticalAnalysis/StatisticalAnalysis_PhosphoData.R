################################################################################
# Statistical analysis of the pYIP samples.
################################################################################

################################################################################
###### Packages
################################################################################
require(reshape2)
require(multcompView) # extract_p()
require(Hmisc) # Pearson correlation
require(matrixStats) # rowMax()

################################################################################
###### Import raw table
################################################################################

BioRep <- c("R1", "R4", "R5")
TimePoints <- c("NS.", "S300.", "S120.", "S30.", "S15.", "S600")
load(file = "RData/05_ParsedTablesNormPSP.RData")

################################################################################
###### Log2-transformation and mean of the technical replicates
################################################################################

################################################################################
# I filter the tables and merge them according to the phosphosite ID (protein accession followed by phospho-position).

ltab2 <- list() ; lmat <- list()
ltab <- lf
for (tab in ltab) {
        tab <- cbind(tab, "psiteID" = tab$phosphoSites)
        
        # Remove non-tyrosines from the pYIP:
        if (grepl("pTyr", names(tab)[12])) {
          tab <- tab[grepl("_Y", tab$psiteID, fixed = T) | grepl("+Y", tab$psiteID, fixed = T),]
        } else {
          tab <- tab[!(grepl("_Y", tab$psiteID, fixed = T) | grepl("+Y", tab$psiteID, fixed = T)),]
        }
        
        mat11 <- tab[,grepl("iRTNorm", names(tab), fixed = T)]
        names(mat11) <- sapply(names(mat11), function(x) strsplit(x, "iRT.", fixed = T)[[1]][2])
        # Remove rows with only missing values:
        n <- sapply(1:nrow(mat11), function(x) length(mat11[x,][is.na(mat11[x,])]))
        tab1 <- tab[n<ncol(mat11),]
        mat11 <- mat11[n<ncol(mat11),]
        row.names(mat11) <- tab1$psiteID
        lmat[[length(lmat)+1]] <- mat11
        ltab2[[length(ltab2)+1]] <- tab1
}
names(ltab2) <- names(lf) ; names(lmat) <- names(lf)

mapping <- sapply(ltab2, function(x) {
  as.character(x$psiteID)
})
mapping <- sort(unique(unlist(mapping)))

df <- data.frame("psiteID" = mapping)
for (i in seq_along(ltab2)) {
  df <- merge(df, ltab2[[i]][,grepl("iRTNorm", names(ltab2[[i]])) | names(ltab2[[i]]) == "psiteID"], by = "psiteID", all = T)
}

keep <- df

################################################################################
# log2-transform the table:
df[,2:ncol(df)] <- log2(df[,2:ncol(df)])
names(df) <- gsub("Intensity.RR_ProFI_", "Log2_", names(df), fixed = T)
names(df) <- gsub("Intensity.ProFI_", "Log2_", names(df), fixed = T)
names(df) <- gsub("Intensity", "Log2_", names(df), fixed = T)
# Correction of naming issues:
names(df) <- gsub("_E.TiO2", "_TiO2", names(df), fixed = T)
names(df) <- gsub("6.TiO2", "_TiO2", names(df), fixed = T)
names(df) <- gsub("E_pTyr", "_pTyr", names(df), fixed = T)
names(df) <- gsub("_EpTyr", "_pTyr", names(df), fixed = T)
names(df) <- gsub("CCF0141_", "CCF01416_", names(df), fixed = T)

keep <- cbind(df, keep)

################################################################################
# For each replicate independently, I calculate the mean of the injection replicates in each time point. I also calculate the standard deviations in order to estimate the noise when replacing missing values.

vec <- names(df)[2:ncol(df)]
vec <- sapply(vec, function(x) {
  strsplit(x, "Inj", fixed = T)[[1]][1]
})
mat <- matrix(ncol = length(unique(vec)), nrow = nrow(df))
colnames(mat) <- unique(vec)
row.names(mat) <- df$psiteID
matsd <- mat
for (i in seq_along(unique(vec))) {
  el <- colnames(mat)[i]
  matel <- df[,grepl(el, names(df), fixed = T)]
  mat[,i] <- rowMeans(matel, na.rm = T)
  matsd[,i] <- sapply(seq_len(nrow(matel)), function(x) {
    sd(matel[x,], na.rm = T)
  })
}
# I calculate the standard deviation for drawing missing values replacement based on the sites under the 10% quantile of the data set:
q5 <- sapply(seq_len(ncol(mat)), function(x) {
  iswithsd <- !is.na(matsd[,x])
  quantile(mat[iswithsd,x], 0.10, na.rm = T)
})
medianStd <- sapply(seq_len(ncol(matsd)), function(x) {
  median(matsd[,x][as.numeric(as.character(mat[,x])) <= q5[x]], na.rm = T)
})
names(medianStd) <- colnames(matsd)

par(mar = c(8, 2, 2, 1))
boxplot(mat, las = 2, main = "Mean of technical repeats")

colnames(mat) <- gsub("Log2", "MeanTech", colnames(mat))
keep <- cbind(keep, mat)

rm(list = c("df", "lf", "lmat", "ltab", "ltab2", "mapping"))
save(keep, file = "RData/07_TableBeforeStat.RData")

################################################################################
###### Statistical analysis
################################################################################

# load(file = "RData/07_TableBeforeStat.RData")
# mat <- keep[,grepl("MeanTech", colnames(keep), fixed = T)]

################################################################################
# For this analysis I generate a loop. This allows the repetition of the entire analysis and we considere regulated the phosphorylation sites that are significantly regulated 90% of the times.
################################################################################

# Create a table that will be used to choose where to replace missing values. The idea is to replace missing values only if there are no values measured in any biological replicate (0 = no replacement, 1 = replacement needed):

BRep <- colnames(mat)
BRep <- gsub("MeanTech_.", "", BRep, fixed = T)
BRep <- gsub("MeanTech_S4_", "", BRep, fixed = T)
BRep <- gsub("MeanTech_S5_", "", BRep, fixed = T)
BRep <- sapply(BRep, function(x) {
  strsplit(x, "_S", fixed = T)[[1]][1]
})
BRep <- sapply(BRep, function(x) {
  strsplit(x, "_N", fixed = T)[[1]][1]
})
BRep <- gsub("CCF[0-9]+", "", BRep)
BRep <- gsub("__", "_", BRep, fixed = T)
UniqueBRep <- sort(unique(BRep))

mat <- mat[,order(BRep)]
BRep <- sort(BRep)

TimePoints <- sapply(colnames(mat), function(x) {
  strsplit(x, "_", fixed = T)[[1]][length(strsplit(x, "_", fixed = T)[[1]])]
})
colnames(mat) <- paste0(BRep, "_", TimePoints)

matRaw <- mat[,order(TimePoints)]

ReplacementMatrix <- matrix(ncol = ncol(matRaw), nrow = nrow(matRaw))
dimnames(ReplacementMatrix) <- dimnames(matRaw)
TimePoints <- sort(TimePoints)

# I replace missing values only the sites with a  minimum of 3 time points with points from a minimum of 2 biological repeats:
#############################################
# For TiO2:
temp <- matRaw[,grepl("TiO2", colnames(matRaw))]
tempnames <- gsub("R[12345]_TiO2_", "", colnames(temp))
tempnames <- gsub("R[12345]_pTyr_", "", tempnames)
temp[!is.na(temp)] <- 1
temp[is.na(temp)] <- 0
colnames(temp) <- tempnames
temp <- sapply(seq_len(nrow(temp)), function(x) {
  table(names(temp[x,][temp[x,]==1])) # Number of non-missing values
})
temp <- lapply(temp, table)
logicalvector <- lapply(temp, function(x) {
  sum(x[names(x) >= 2]) >= 3
}) 
logicalvector <- unlist(logicalvector)
ReplacementAction <- (sapply(temp, length) > 0) & logicalvector
names(ReplacementAction) <- row.names(matRaw)
ReplacementActionTiO2 <- ReplacementAction

# For pYIP:
temp <- matRaw[,grepl("pTyr", colnames(matRaw))]
tempnames <- gsub("R[12345]_TiO2_", "", colnames(temp))
tempnames <- gsub("R[12345]_pTyr_", "", tempnames)
temp[!is.na(temp)] <- 1
temp[is.na(temp)] <- 0
colnames(temp) <- tempnames
temp <- sapply(seq_len(nrow(temp)), function(x) {
  table(names(temp[x,][temp[x,]==1])) # Number of non-missing values
})
temp <- lapply(temp, table)
logicalvector <- lapply(temp, function(x) {
  sum(x[names(x) >= 2]) >= 3
}) 
logicalvector <- unlist(logicalvector)
ReplacementAction <- (sapply(temp, length) > 0) & logicalvector
names(ReplacementAction) <- row.names(matRaw)
ReplacementActionpTyr <- ReplacementAction


vec <- gsub("R[1-5]", "", colnames(matRaw))
for (el in unique(vec)) {
  matel <- matRaw[,vec == el]
  sumNA <- sapply(seq_len(nrow(matel)), function(x) {
    length(matel[x,][is.na(matel[x,])])
  })
  if (grepl("pTyr", el)) {
    ReplacementMatrix[sumNA == ncol(matel) & ReplacementActionpTyr, vec == el] <- 1
  } else {
    ReplacementMatrix[sumNA == ncol(matel) & ReplacementActionTiO2, vec == el] <- 1
  }
}
ReplacementMatrix[is.na(ReplacementMatrix)] <- 0
rm(ReplacementAction)

# Will only loop where replacement is needed:
TFloop <- rowSums(ReplacementMatrix) > 0

# For replacement of missing values:
## I replace missing values with a random value generated around the 5% quantile of the intensities for each condition, using the standard deviation of the technical repeats.
noise <- sapply(seq_len(ncol(matRaw)), function(x) {
  quantile(matRaw[,x], 0.05, na.rm = T)
})
# Preparing the standard deviations for the replacement:
sdnames <- names(medianStd)
sdnames <- gsub("Log2_.", "", sdnames, fixed = T)
sdnames <- gsub("Log2_S4_", "", sdnames, fixed = T)
sdnames <- gsub("Log2_S5_", "", sdnames, fixed = T)
sdnames <- gsub("CCF[0-9]+", "", sdnames)
sdnames <- gsub("__", "_", sdnames, fixed = T)
names(medianStd) <- sdnames
medianStd <- medianStd[match(colnames(matRaw), names(medianStd))]

################################################################################

nLoops <- 200
MatBeforeStat <- vector(mode = "list")
MatAfterStat <- vector(mode = "list")
while (nLoops > 0) {
  # keep only the rows where values need to be replaced:
  if (nLoops == 1) {
    matTF <- matRaw
  } else {
    matTF <- matRaw[TFloop,]
  }
  mat2 <- matTF
  # REPLACEMENT VALUE:
  for (i in seq_len(ncol(mat2))) {
    vrep <- rnorm(n = length(matTF[,i]),mean = noise[i], sd = as.numeric(as.character(medianStd[i])))
    if (nLoops == 1) {
      mat2[,i][ReplacementMatrix[,i]==1] <- vrep[ReplacementMatrix[,i]==1]
    } else {
      mat2[,i][ReplacementMatrix[TFloop,i]==1] <- vrep[ReplacementMatrix[TFloop,i]==1]
    }
  }
  
  if (nLoops == 1) {
    boxplot(matTF, main = "Raw values", las = 2)
    boxplot(mat2, main = "After replacement of missing values", las = 2)
  }
  
  # For each phosphorylation site, I calculate the mean value per biological replicate after replacement of missing values. The mean of these values are used for normalisation.
  RowNorm <- rowMeans(mat2, na.rm = T)
  
  vec <- UniqueBRep
  matmedkin <- matrix(ncol = length(unique(vec)), nrow = nrow(mat2))
  for (i in seq_along(unique(vec))) {
    matel <- mat2[,grepl(unique(vec)[i], colnames(mat2), fixed = T)]
    matmedkin[,i] <- sapply(seq_len(nrow(matel)), function(x) {
      mean(as.numeric(as.character(matel[x,])), na.rm = T)
    })
  }
  colnames(matmedkin) <- unique(vec)
  MeanKinetics <- sapply(seq_len(nrow(matmedkin)), function(x) {
    mean(matmedkin[x,], na.rm = T)
  })
  # Calculation of the normalisation factors:
  for (i in seq_len(nrow(matmedkin))) {
    matmedkin[i,] <- matmedkin[i,] / MeanKinetics[i]
  }
  
  for (i in seq_along(unique(vec))) {
    matel <- mat2[,grepl(unique(vec)[i], colnames(mat2))]
    for (j in seq_len(nrow(matel))) {
      matel[j,] <- matel[j,] / matmedkin[j,colnames(matmedkin) == unique(vec)[i]] 
    }
    mat2[,grepl(unique(vec)[i], colnames(mat2))] <- matel
  }
  if (nLoops == 1) {
    boxplot(mat2, las = 2, main = "After replicates re-alignment")
  }

  ########################################################################
  # I perform an anova on each phosphorylation site across the 6 time points (only the phosphorylation sites with a minimum of three quantification values in a minimum of two biological replicates, the same time points being considered). Then, I correct the pvalue with Tukey HSD.
  
  df <- mat2
  df <- melt(df)
  TimePoints <- sort(unique(TimePoints))
  vec <- ifelse(grepl(TimePoints[1], df$Var2, fixed = T), "T0", NA)
  vec[grepl(TimePoints[2], df$Var2, fixed = T)] <- "T15"
  vec[grepl(TimePoints[3], df$Var2, fixed = T)] <- "T30"
  vec[grepl(TimePoints[4], df$Var2, fixed = T)] <- "T120"
  vec[grepl(TimePoints[5], df$Var2, fixed = T)] <- "T300"
  vec[grepl(TimePoints[6], df$Var2, fixed = T)] <- "T600"
  df <- data.frame(df, "TimePoint" = vec)
  
  df <- df[!is.na(df$value),]
  
  # For calculation of the mean of all loops:
  MatBeforeStat[[length(MatBeforeStat)+1]] <- mat2[,order(colnames(mat2))]
  
  psites <- as.character(unique(row.names(mat2)))
  
  anov <- vector(mode = "list")
  npTuk <- sort(c("T15-T120", "T30-T120", "T300-T120", "T600-T120", "T30-T15", "T300-T15", "T600-T15", "T300-T30", "T600-T30", "T600-T300", "T15-T0", "T30-T0", "T120-T0", "T300-T0", "T600-T0"))
  
  for (i in 1:length(psites)) {
    el <- psites[i]
    sub <- subset(df, df$Var1==el)
    subtemp <- sub[!is.na(sub$value),]
    # if (length(table(subtemp$TimePoint)[table(sub$TimePoint) >= 2]) >= 1 & length(unique(sub$TimePoint)) > 3) {
    # I keep only the sites with a  minimum of 3 time points with points from a minimum of 2 biological repeats:
    #############################################
    temp <- table(subtemp$TimePoint)
    temp <- table(temp)
    if (sum(temp[names(temp) >= 2]) >= 3) {
      anovA <- aov(sub$value~sub$TimePoint, data = sub)
      anovB <- anova(anovA)
      pAnova <- anovB$"Pr(>F)"[1]
      Tuk <- TukeyHSD(anovA)
      pTuk <- extract_p(Tuk)[order(names(extract_p(Tuk)))]
      if (length(pTuk[[1]]) == 15) {
        vec <- c("psiteID" = el, "pAnova"=pAnova, pTuk[[1]][match(npTuk, names(pTuk[[1]]))])
        anov[[length(anov)+1]] <- vec
      } else {
        pTvec <- pTuk[[1]]
        pTvec <- c(pTvec, rep(NA, 15-length(pTvec)))
        names(pTvec)[(length(pTuk[[1]]) + 1):length(pTvec)] <- setdiff(npTuk, names(pTvec))
        vec <- c("psiteID" = el, "pAnova"=pAnova, pTvec[match(npTuk, names(pTvec))])
        anov[[length(anov)+1]] <- vec
      }
    } 
  } 
    
  anov <- lapply(anov, function(x) x[names(anov[[1]])])
  mat <- rbind(anov[[1]], anov[[2]])
  for (i in 3:length(anov)) {
    mat <- rbind(mat, anov[[i]])
  }
  row.names(mat) <- mat[,1]
  matk <- mat
  
  
  # I calculate the mean values of the biological replicates. I also keep track of the minimum number of point per comparison to discard the values with less than 2 points in one of the two conditions:
  mat <- matrix(ncol = 15, nrow = nrow(mat2))
  matnum <- matrix(ncol = 15, nrow = nrow(mat2))
  dimnames(mat)[[2]] <- c( "T120-T0", "T600-T120", "T300-T0", "T300-T15", "T15-T0", "T30-T0", "T300-T30", "T600-T0", "T15-T120", "T600-T15", "T600-T300", "T600-T30", "T30-T120", "T300-T120", "T30-T15")
  dimnames(matnum)[[2]] <- c( "T120-T0", "T600-T120", "T300-T0", "T300-T15", "T15-T0", "T30-T0", "T300-T30", "T600-T0", "T15-T120", "T600-T15", "T600-T300", "T600-T30", "T30-T120", "T300-T120", "T30-T15")
  for (i in 1:nrow(mat2)) {
    sub <- melt(mat2[i,], na.rm = T)
    mat[i,1] <- mean(sub[grepl("120", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("NS.", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,1] <- min(length(sub[grepl("120", rownames(sub), fixed = T),1]), length(sub[grepl("NS.", rownames(sub), fixed = T),1]))
    mat[i,2] <- mean(sub[grepl("600", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("120", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,2] <- min(length(sub[grepl("600", rownames(sub), fixed = T),1]), length(sub[grepl("120", rownames(sub), fixed = T),1]))
    mat[i,3] <- mean(sub[grepl("300", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("NS.", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,3] <- min(length(sub[grepl("300", rownames(sub), fixed = T),1]), length(sub[grepl("NS.", rownames(sub), fixed = T),1]))
    mat[i,4] <- mean(sub[grepl("300", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("15.", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,4] <- min(length(sub[grepl("300", rownames(sub), fixed = T),1]), length(sub[grepl("15.", rownames(sub), fixed = T),1]))
    mat[i,5] <- mean(sub[grepl("15.", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("NS.", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,5] <- min(length(sub[grepl("15.", rownames(sub), fixed = T),1]), length(sub[grepl("NS.", rownames(sub), fixed = T),1]))
    mat[i,6] <- mean(sub[grepl("30.", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("NS.", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,6] <- min(length(sub[grepl("30.", rownames(sub), fixed = T),1]), length(sub[grepl("NS.", rownames(sub), fixed = T),1]))
    mat[i,7] <- mean(sub[grepl("300", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("30.", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,7] <- min(length(sub[grepl("300", rownames(sub), fixed = T),1]), length(sub[grepl("30.", rownames(sub), fixed = T),1]))
    mat[i,8] <- mean(sub[grepl("600", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("NS.", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,8] <- min(length(sub[grepl("600", rownames(sub), fixed = T),1]), length(sub[grepl("NS.", rownames(sub), fixed = T),1]))
    mat[i,9] <- mean(sub[grepl("15.", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("120", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,9] <- min(length(sub[grepl("15.", rownames(sub), fixed = T),1]), length(sub[grepl("120", rownames(sub), fixed = T),1]))
    mat[i,10] <- mean(sub[grepl("600", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("15.", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,10] <- min(length(sub[grepl("600", rownames(sub), fixed = T),1]), length(sub[grepl("15.", rownames(sub), fixed = T),1]))
    mat[i,11] <- mean(sub[grepl("600", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("300", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,11] <- min(length(sub[grepl("600", rownames(sub), fixed = T),1]), length(sub[grepl("300", rownames(sub), fixed = T),1]))
    mat[i,12] <- mean(sub[grepl("600", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("30.", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,12] <- min(length(sub[grepl("600", rownames(sub), fixed = T),1]), length(sub[grepl("30.", rownames(sub), fixed = T),1]))
    mat[i,13] <- mean(sub[grepl("30.", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("120", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,13] <- min(length(sub[grepl("30.", rownames(sub), fixed = T),1]), length(sub[grepl("120", rownames(sub), fixed = T),1]))
    mat[i,14] <- mean(sub[grepl("300", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("120", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,14] <- min(length(sub[grepl("300", rownames(sub), fixed = T),1]), length(sub[grepl("120", rownames(sub), fixed = T),1]))
    mat[i,15] <- mean(sub[grepl("30.", rownames(sub), fixed = T),], na.rm = T) - mean(sub[grepl("15.", rownames(sub), fixed = T),], na.rm = T)
    matnum[i,15] <- min(length(sub[grepl("30.", rownames(sub), fixed = T),1]), length(sub[grepl("15.", rownames(sub), fixed = T),1]))
  }
  # Remove FC when there is only one repeat:
  mat[matnum <2] <- NA
  mat <- mat[,order(colnames(mat))]
  row.names(mat) <- row.names(mat2)
  
  # Retrieve the best FC:
  matFC <- mat
  matpT <- apply(matk[match(row.names(mat), matk[,1]),3:ncol(matk)], 2, as.numeric)
  matFC <- matFC[,order(colnames(matFC))]
  matpT <- matpT[,order(colnames(matpT))]
  k <- matpT<=0.05
  k[is.na(k)] <- FALSE
  k[is.na(matFC)] <- FALSE
  matFC[!k] <- NA
  btp <- sapply(1:nrow(matFC), function(x) as.numeric(matFC[x,][order(abs(matFC[x,]), decreasing = T)][1]))
  BestTimePoint <- sapply(1:nrow(matFC), function(x) names(matFC[x,][order(abs(matFC[x,]), decreasing = T)][1]))
  BestTimePoint[is.na(btp)] <- NA
  
  colnames(mat) <- paste0(colnames(mat), ".FC")
  mat <- cbind(mat, "BestTimePoint"=BestTimePoint, "BestFC"=btp)
  df <- data.frame("psiteID" = row.names(mat2), mat[match(row.names(mat2), row.names(mat)),])
  
  mer <- merge(df, matk, by = "psiteID", all = T)
  
  
  # SIGNIFICATIVITY THRESHOLD:
  mer$Regulation <- ifelse(!is.na(mer$BestFC) & abs(as.numeric(as.character(mer$BestFC))) > log2(1.5) & (as.numeric(as.character(mer$pAnova))<=0.05), "Regulated", "Not regulated")
  
  MatAfterStat[[length(MatAfterStat)+1]] <- mer
  
  nLoops <- (nLoops - 1)
  print(nLoops)
}

################################################################################
# Merging of the result of the loops:

nLoops <- length(MatAfterStat)

# Phosphosites significantly regulated across 90% of the statistical tests performed:
reg <- ifelse(MatAfterStat[[nLoops]]$Regulation == "Regulated", 1, 0)
names(reg) <-  MatAfterStat[[nLoops]]$psiteID 
regi <- reg
for (i in seq_len(nLoops-1)) {
  vec <- ifelse(MatAfterStat[[i]]$Regulation == "Regulated", 1, 0)[match(names(regi), MatAfterStat[[i]]$psiteID)]
  reg[!is.na(vec) & !is.na(reg)] <- reg[!is.na(vec) & !is.na(reg)] + vec[!is.na(vec) & !is.na(reg)]
  reg[is.na(vec) & regi == 1 & !is.na(reg)] <- as.numeric(reg[is.na(vec) & regi == 1 & !is.na(reg)]) + 1
}
threshold <- .9*nLoops # 90% of the draws
regval <- reg
reg <- regval>=threshold
regval <- regval/nLoops

matmean <- MatBeforeStat[[nLoops]]
for (i in seq_len(nLoops-1)) {
  mati <- MatBeforeStat[[i]][match(row.names(matmean), row.names(MatBeforeStat[[i]])),]
  for (j in seq_len(ncol(matmean))) {
    matmean[,j] <- sapply(seq_len(nrow(matmean)), function(x) {
      mean(c(matmean[x,j], mati[x,j]), na.rm = T)
    })
  }
}

colnames(matmean) <- paste0("MeanLoops_", colnames(matmean))
keep <- cbind(keep, matmean)

matmean <- apply(MatAfterStat[[nLoops]][,c(2:16, 19:34)], 2, as.character)
matmean <- apply(matmean, 2, as.numeric)
row.names(matmean) <- MatAfterStat[[nLoops]]$psiteID
for (i in seq_len(nLoops-1)) {
  mati <- MatAfterStat[[i]][match(row.names(matmean), as.character(MatAfterStat[[i]]$psiteID)),]
  mati <- apply(mati, 2, as.character)[,c(2:16, 19:34)]
  mati <- apply(mati, 2, as.numeric)
  for (j in seq_len(ncol(matmean))) {
    matmean[,j] <- sapply(seq_len(nrow(matmean)), function(x) {
      mean(c(matmean[x,j], mati[x,j]), na.rm = T)
    })
  }
}

# Retrieve the best FC:
matFC <- abs(matmean[,grepl(".FC", colnames(matmean), fixed = T)])
matpT <- matmean[,17:31]
matFC <- matFC[,order(colnames(matFC))]
matpT <- matpT[,order(colnames(matpT))]
colnames(matFC) <- gsub(".", "-", colnames(matFC), fixed = T)
matpT <- matpT[,match(gsub("-FC", "", colnames(matFC)), colnames(matpT))]
k <- matpT<=0.05
k[is.na(k)] <- FALSE
matFC[!k] <- NA
btp <- sapply(1:nrow(matFC), function(x) as.numeric(matFC[x,][order(abs(matFC[x,]), decreasing = T)][1]))
BestTimePoint <- sapply(1:nrow(matFC), function(x) names(matFC[x,][order(abs(matFC[x,]), decreasing = T)][1]))
BestTimePoint[is.na(btp)] <- NA

finalMat <- data.frame(matmean, "BestTimePoint"=BestTimePoint, "BestFC"=btp, "ProportioPassStat" = regval[match(row.names(matmean), names(regval))], "Regulation" = reg[match(row.names(matmean), names(reg))], "psiteID" = row.names(matmean))
finalMat$Regulation[is.na(finalMat$Regulation)] <- "Not tested"

keep <- keep[,c(1:135, 137:ncol(keep))]
export <- merge(keep, finalMat, by = "psiteID", all = T)

load(file = "RData/03_psiteInfoPhosphoSitePlus.RData")

# psiteInfo[names(psiteInfo) == "Q9Z0R6_Y922"][[1]]$GeneID # Itsn2
vec <- vector(mode = "character")
for (i in seq_len(nrow(export))) {
  if (length(unique(psiteInfo[names(psiteInfo) == export$psiteID[i]][[1]]$GeneID))==1) {
    vec[i] <- unique(psiteInfo[names(psiteInfo) == export$psiteID[i]][[1]]$GeneID)
  } else {
    paste(unique(psiteInfo[names(psiteInfo) == export$psiteID[i]][[1]]$GeneID), collapse = "|")
  }
}
export$GeneID  <- vec

write.table(export, file = "StatisticalAnalysis/PhosphoStatResults.txt", row.names = F, sep = "\t", quote = F)

RegulatedHits <- as.character(export$psiteID[export$Regulation == "TRUE"])
RegulatedHits <- RegulatedHits[!is.na(RegulatedHits)]
save(RegulatedHits, file = "RData/RegulatedpSites.Rdata")
save(export, file = "RData/08_StatResults_Phospho.RData")

# load(file = "RData/RegulatedpSites.Rdata")
# load(file = "RData/08_StatResults_Phospho.RData")

sessionInfo()
