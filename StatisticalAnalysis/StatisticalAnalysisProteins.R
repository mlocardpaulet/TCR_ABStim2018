################################################################################
# Correlation analysis of the protein samples.
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

BioRep <- c("R1", "R3", "R4", "R5")
TimePoints <- c("NS.", "S300.", "S120.", "S30.", "S15.", "S600")
load(file = "RData/06_Prottab.RData")

################################################################################
# Filter the table:

prot <- prot[!grepl("CON_", prot$Protein.IDs, fixed = T)& !grepl("REV_", prot$Protein.IDs, fixed = T),]

################################################################################
###### Log2-transformation and mean of the technical replicates
################################################################################

################################################################################
################################################################################
# For each replicate independently, I calculate the mean of the injection replicates in each time point:


lmat2 <- list()
for (el1 in BioRep) {
        tab <- prot[,grepl(paste0("_", el1, "_"), names(prot), fixed = T) & grepl("LFQ.intensity.", names(prot), fixed = T)]
        tab[tab==0] <- NA
        tab <- cbind(prot$Majority.protein.IDs, log2(tab))
        lmat2[[length(lmat2)+1]] <- tab
}
names(lmat2) <- BioRep

lm <- list()
lmsd <- list()
for (i in 1:length(lmat2)) {
        tab <- lmat2[[i]]
        mat <- matrix(ncol = length(TimePoints), nrow = nrow(tab))
        dimnames(mat)[[2]] <- TimePoints
        dimnames(mat)[[1]] <- tab[,1]
        matsd <- mat
        j = 1
        for (el in TimePoints) {
                sub <- tab[,grepl(el, names(tab), fixed = T)]
                mat[,j] <- rowMeans(sub, na.rm = T)
                matsd[,j] <- sapply(seq_len(nrow(sub)), function(x) {
                  sd(sub[x,], na.rm = T)
                })
                j <- j+1
        }
        lm[[length(lm)+1]] <- mat
        lmsd[[length(lmsd)+1]] <- matsd
}
names(lm) <- BioRep
names(lmsd) <- BioRep

lmsd2 <- vector(mode ="list")
for (i in seq_along(lmsd)) {
  matsd <- lmsd[[i]]
  mat <- lm[[i]]
  q5 <- sapply(seq_len(ncol(mat)), function(x) {
    iswithsd <- !is.na(matsd[,x])
    quantile(mat[iswithsd,x], 0.10, na.rm = T)
  })
  medianStd <- sapply(seq_len(ncol(matsd)), function(x) {
    median(matsd[,x][as.numeric(as.character(mat[,x])) <= q5[x]], na.rm = T)
  })
  names(medianStd) <- colnames(matsd)
  lmsd2[[i]] <- medianStd
}

names(lmsd2) <- names(lm)
rm(lmsd)

lm <- lapply(lm, function(x) x[,c(1, 5, 4, 3, 2, 6)])
lmsd2 <- lapply(lmsd2, function(x) x[c(1, 5, 4, 3, 2, 6)])

# Will only loop where replacement is needed:
TFloop <- vector(mode = "list")
for (i in seq_along(lm)) {
  mat <- lm[[i]]
  TFloop[[i]] <- sapply(seq_len(nrow(mat)), function(x) {
    length(mat[x,][is.na(mat[x,])]) > 0
  })
}

################################################################################
###### Statistical analysis
################################################################################

################################################################################
# For this analysis I generate a loop. This allows the repetition of the entire analysis and we considere regulated the phosphorylation sites that are significantly regulated 90% of the times.
################################################################################


nLoops <- 200
finalListAnova <- list()
set.seed(231)
while (nLoops > 0) {
        
        # I replace missing values with a random value generated around the 5% quantile of the intensities for each condition, using the standard deviation of the corresponding condition (standard deviation generated with the values included in the 5% quantile only). I don't need any normalisation since I work with LFQ values.
        
        lmat2 <- list()
        lmat3 <- list()
        par(mfrow=c(1, length(lm)))
        for (j in 1:length(lm)) {
          mat <- lm[[j]]
          noisevec <- sapply(seq_len(ncol(mat)), function(x) {
            quantile(mat[,x], 0.05, na.rm = T)
          })
          if (nLoops > 1) {
            mat <- mat[TFloop[[j]],]
          }
          stdv <- lmsd2[[j]]
          mat2 <- mat
          for (i in seq_len(ncol(mat))) {
            noise <- noisevec[i]
            stdvi <- stdv[i]
            vrep <- rnorm(length(mat[,i]), mean = noise, sd = stdvi)
            mat2[,i][is.na(mat[,i])] <- vrep[is.na(mat[,i])]
          }
          lmat2[[length(lmat2)+1]] <- mat2
          val <- rowMeans(mat2, na.rm = T)
          lmat3[[length(lmat3)+1]] <- data.frame(mat, "NormalisationValue"=val)
          if (nLoops == 1) {
            title <- paste0(names(lm)[j], ": quan values after replacement MV")
            boxplot(mat2, main = title, las = 2)
          }
        }
        names(lmat2) <- names(lm) # lmat2 contains the tables with no missing values.
        names(lmat3) <- names(lm) # lmat3 contains the raw values and normalisation values.
        
        # For each protein, I calculate the mean value per biological replicate after replacement of missing values. The mean of these values are used for normalisation.
        
        l <- list()
        for (i in 1:length(lmat3)) {
                vec <- lmat3[[i]][,7]
                na <- dimnames(lmat3[[i]])[[1]]
                l[[length(l)+1]] <- cbind(na, vec)
        }
        names(l) <- names(lm)
        normtab <- rbind(l[[1]], l[[2]])
        for (i in 3:length(l)) {
                normtab <- rbind(normtab, l[[i]])
        }
        normtab <- as.data.frame(normtab, stringsAsFactors = F)
        normtab[,2] <- as.numeric(normtab[,2])
        NormVal1 <- lapply(unique(normtab$na), function(x) normtab[normtab$na==x,2])
        names(NormVal1) <- unique(normtab$na)
        NormVal2 <- sapply(NormVal1, function(x) {
          mean(as.numeric(as.character(x)), na.rm = T)
        })
        l2 <- l
        par(mfrow=c(1,length(l)))
        for (i in 1:length(l)) {
                tab <- as.data.frame(l[[i]], stringsAsFactors = F)
                normval <- NormVal2[match(tab$na, names(NormVal1))]
                vec2 <- as.numeric(tab$vec)-normval
                l2[[i]] <- data.frame(tab, "NormValue"=vec2)
                #boxplot(vec2, main = paste0("Normalisation values of ", names(l)[i]))
        }
        
        lmat4 <- lmat2
        for (i in 1:length(lmat4)) {
                mat <- lmat4[[i]][,1:6]
                matn <- l2[[i]]
                for (j in 1:nrow(mat)) {
                        mat[j,] <- mat[j,]-matn[j,3]
                }
                lmat3[[i]] <- mat
                #boxplot(mat, main = paste0("Normalised values of ", names(l)[i]))
        }
        
        
        
        
        # I combine the tables:

        lmat4 <- lmat3
        for (i in 1:length(lmat3)) {
                tab <- lmat3[[i]][,1:6]
                dimnames(tab)[[2]] <- paste0(names(lm)[i], dimnames(tab)[[2]])
                # tab <- data.frame("Majority.protein.IDs" = dimnames(lm[[i]])[[1]], tab)
                tab <- data.frame("Majority.protein.IDs" = row.names(tab), tab)
                lmat4[[i]] <- tab
        }
        names(lmat4) <- names(lmat3)
        mer <- merge(prot, lmat4[[1]], by = "Majority.protein.IDs", all = T)
        for (i in 2:length(lmat4)) {
                mer <- merge(mer, lmat4[[i]], by = "Majority.protein.IDs", all = T)
        }
        # I remove the rows with less than 12 values in the columns 684:707:
        k <- sapply(1:nrow(mer), function(x) {length(mer[x,684:707][!is.na(mer[x,684:707])])})
        mer <- mer[k > 12,]
        
        
        ########################################################################
        # I perform an anova on each phosphorylation site across the 6 time points. Then, I correct the pvalue with Tukey HSD.
        
        df <- data.frame(mer$Majority.protein.IDs, mer[,684:707])
        df <- melt(df)
        vec <- ifelse(grepl(TimePoints[1], df$variable, fixed = T), "T0", NA)
        vec[grepl(TimePoints[2], df$variable, fixed = T)] <- "T300"
        vec[grepl(TimePoints[3], df$variable, fixed = T)] <- "T120"
        vec[grepl(TimePoints[4], df$variable, fixed = T)] <- "T30"
        vec[grepl(TimePoints[5], df$variable, fixed = T)] <- "T15"
        vec[grepl(TimePoints[6], df$variable, fixed = T)] <- "T600"
        df <- data.frame(df, "TimePoint" = vec)
        
        psites <- as.character(unique(df$mer.Majority.protein.IDs))
        
        anov <- list()
        
        for (i in 1:length(psites)) {
                el <- psites[i]
                sub <- subset(df, df$mer.Majority.protein.IDs==el)
                anovA <- aov(sub$value~sub$TimePoint, data = sub)
                anovB <- anova(anovA)
                pAnova <- anovB$"Pr(>F)"[1]
                Tuk <- TukeyHSD(anovA)
                pTuk <- extract_p(Tuk)
                BestpTuk <- min(as.numeric(pTuk[[1]])) # for a best ptuk corresponding to the minimal pTukey
                btp <- abs(Tuk$`sub$TimePoint`[,1])[Tuk$`sub$TimePoint`[,4]<=0.05]
                if (length(btp)==0) { 
                        btp <- NA 
                        BestTimePoint <- NA
                } else {
                        btp <- sort(btp, decreasing = T)[1]
                        BestTimePoint <- names(btp)
                }
                vec <- c("pAnova"=pAnova, pTuk[[1]], "BestpTuk" = BestpTuk, "BestTimePoint" = BestTimePoint, "BestFC" = as.numeric(btp))
                anov[[length(anov)+1]] <- vec
        }
        names(anov) <- psites
        
        anov <- lapply(anov, function(x) x[names(anov[[1]])])
        mat <- rbind(anov[[1]], anov[[2]])
        for (i in 3:length(anov)) {
                mat <- rbind(mat, anov[[i]])
        }
        df <- data.frame("Majority.protein.IDs"=psites, mat)
        mer2 <- merge(mer[,c(1,684:707)], df, by = "Majority.protein.IDs", all = T)
        
        # I calculate the mean values of the biological replicates:
        mat <- matrix(ncol = 15, nrow = nrow(mer2))
        dimnames(mat)[[2]] <- c( "T120-T0", "T600-T120", "T300-T0", "T300-T15", "T15-T0", "T30-T0", "T300-T30", "T600-T0", "T15-T120", "T600-T15", "T600-T300", "T600-T30", "T30-T120", "T300-T120", "T30-T15")
        for (i in 1:nrow(mer2)) {
                sub <- mer2[i,2:25]
                mat[i,1] <- mean(as.numeric(sub[,grepl("120", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("NS.", names(sub), fixed = T)]), na.rm = T)
                mat[i,2] <- mean(as.numeric(sub[,grepl("600", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("120", names(sub), fixed = T)]), na.rm = T)
                mat[i,3] <- mean(as.numeric(sub[,grepl("300", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("NS.", names(sub), fixed = T)]), na.rm = T)
                mat[i,4] <- mean(as.numeric(sub[,grepl("300", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("15.", names(sub), fixed = T)]), na.rm = T)
                mat[i,5] <- mean(as.numeric(sub[,grepl("15.", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("NS.", names(sub), fixed = T)]), na.rm = T)
                mat[i,6] <- mean(as.numeric(sub[,grepl("30.", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("NS.", names(sub), fixed = T)]), na.rm = T)
                mat[i,7] <- mean(as.numeric(sub[,grepl("300", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("30.", names(sub), fixed = T)]), na.rm = T)
                mat[i,8] <- mean(as.numeric(sub[,grepl("600", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("NS.", names(sub), fixed = T)]), na.rm = T)
                mat[i,9] <-  mean(as.numeric(sub[,grepl("15.", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("120", names(sub), fixed = T)]), na.rm = T)
                mat[i,10] <- mean(as.numeric(sub[,grepl("600", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("15.", names(sub), fixed = T)]), na.rm = T)
                mat[i,11] <- mean(as.numeric(sub[,grepl("600", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("300", names(sub), fixed = T)]), na.rm = T)
                mat[i,12] <- mean(as.numeric(sub[,grepl("600", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("30.", names(sub), fixed = T)]), na.rm = T)
                mat[i,13] <- mean(as.numeric(sub[,grepl("30.", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("120", names(sub), fixed = T)]), na.rm = T)
                mat[i,14] <- mean(as.numeric(sub[,grepl("300", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("120", names(sub), fixed = T)]), na.rm = T)
                mat[i,15] <- mean(as.numeric(sub[,grepl("30.", names(sub), fixed = T)]), na.rm = T) - mean(as.numeric(sub[,grepl("15.", names(sub), fixed = T)]), na.rm = T)
        }
        
        df <- data.frame("Majority.protein.IDs" = mer2$Majority.protein.IDs, mat)
        names(df)[2:ncol(df)] <- paste0(names(df)[2:ncol(df)], ".FC")
        mer2 <- merge(mer2, df, by = "Majority.protein.IDs", all = T)
        
        BestFC <- rep(NA, nrow(mer2))
        for (i in 1:nrow(mer2)) {
                tp <- as.character(mer2$BestTimePoint[i])
                if (!(is.na(tp))) {
                        tp <- sub("-", ".", tp)
                        BestFC[i] <- mer2[i,grepl(paste0(tp, ".FC"), names(mer2))]
                }
        }
        mer2 <- cbind(mer2, "BestFC"=BestFC)
        # SIGNIFICATIVITY THRESHOLD:
        mer2$Regulation <- ifelse(!is.na(mer2$BestpTuk) & (as.numeric(as.character(mer2$BestpTuk))<=0.05) & abs(as.numeric(as.character(mer2$BestFC))>=log2(1.75) & (as.numeric(as.character(mer2$pAnova))<=0.05)), "Regulated", "Not regulated")
        
        finalListAnova[[length(finalListAnova)+1]] <- mer2
        
        nLoops <- (nLoops - 1)
        print(nLoops)
}

################################################################################
# Merging of the result of the loops:

nLoops <- length(finalListAnova)

# Phosphosites significantly regulated across 90% of the statistical tests performed:
reg <- lapply(finalListAnova, function(x) x[,"Regulation"])
for (i in seq_along(reg)) {
  names(reg[[i]]) <- finalListAnova[[i]]$Majority.protein.IDs
}
threshold <- .9*nLoops # 90% of the draws
na <- finalListAnova[[length(finalListAnova)]]$Majority.protein.IDs
reg1 <- data.frame("Majority.protein.IDs" = na, "1" = ifelse(reg[[1]][match(na, names(reg[[1]]))] == "Regulated", 1, 0))
for (i in (2 : length(reg))) {
  reg1 <- cbind(reg1, ifelse(reg[[i]][match(na, names(reg[[i]]))] == "Regulated", 1, 0))
}
names(reg1)[3:ncol(reg1)] <-  as.numeric(2 : length(reg))
# reg1[is.na(reg1)] <- 0
reg <- reg1
reg1 <- rowSums(reg[2:ncol(reg)])
reg1[is.na(reg$X1)] <- reg[is.na(reg$X1), (nLoops + 1)] * nLoops
reg <- reg1>=threshold & !is.na(reg1)
reg1 <- reg1/nLoops
names(reg1) <- names(reg)


lmat <- lapply(finalListAnova, function(x) x[,c(2:42, 44:59)])
lmat <- lapply(lmat, apply, 2, as.numeric)
for (i in seq_along(lmat)) {
  row.names(lmat[[i]]) <- finalListAnova[[i]]$Majority.protein.IDs
}
finalMat <- matrix(ncol = ncol(lmat[[nLoops]]), nrow = nrow(lmat[[nLoops]]))
colnames(finalMat) <- colnames(lmat[[1]])
row.names(finalMat) <- finalListAnova[[nLoops]]$Majority.protein.IDs
source("RFunctions/RBindList.R")
for (i in 1:nrow(finalMat)) {
  proti <- row.names(finalMat)[i]
  mat <- lapply(lmat, function(x) {
    x[row.names(x) == proti,]
  })
  mat <- t(RBindList(mat))
  finalMat[i,] <- rowMeans(mat, na.rm = T)
}
finalMat <- data.frame("Majority.protein.IDs" = row.names(finalMat), finalMat, "Regulation" = reg[match(row.names(finalMat), names(reg))])

finalMat <- data.frame(finalMat, "FrequenceSignificance" = as.numeric(reg1)[match(as.character(finalMat$Majority.protein.IDs), names(reg1))])
finalMat <- finalMat[,!grepl("BestpTuk", names(finalMat), fixed = T) & !grepl("BestFC", names(finalMat), fixed = T)]

# Retrieve the best FC:
matFC <- abs(finalMat[,grepl(".FC", names(finalMat), fixed = T)])
matpT <- finalMat[,27:41]
matpT <- matpT[,match(gsub(".FC", "", names(matFC)), names(matpT))]
k <- matpT<=0.05
k[is.na(k)] <- FALSE
matFC[!k] <- NA
btp <- sapply(1:nrow(matFC), function(x) sort(as.numeric(matFC[x,]), decreasing = T)[1])
BestTimePoint <- sapply(1:nrow(matFC), function(x) names(sort(matFC[x,], decreasing = T))[1])
val <- rowMins(as.matrix(matpT), na.rm = T)
finalMat <- data.frame(finalMat, "BestTimePoint"=BestTimePoint, "BestFC"=btp, "BestpTuk"=val)



export <- merge(prot, finalMat, by = "Majority.protein.IDs", all = T)

#export <- data.frame(cbind(export, "MaxLocalisationScore" = mer$MaxLocalisationScore))

# Issues with openning the table in excel. There are some "\n":
# mat <- export[,673:683]
# mat <- apply(mat, 2, as.character)
# mat <- gsub("\n", "", mat, fixed = T)
# mat <- gsub("\t", "", mat, fixed = T)
# mat <- gsub(";", "|", mat, fixed = T)
# mat <- gsub(",", "|", mat, fixed = T)
# export[,673:683] <- mat
export <- export[,c(1:672,684:743)]
write.table(export, "SupTables/ProtfinalStat.txt", row.names = F, sep = "\t", quote = F)
write.csv(export, "SupTables/ProtfinalStat.csv", row.names = F)

RegulatedHits <- as.character(export$Majority.protein.IDs[export$Regulation==TRUE])
RegulatedHits <- RegulatedHits[!is.na(RegulatedHits)]
save(RegulatedHits, file = "RData/RegulatedProt.Rdata")

sessionInfo()
