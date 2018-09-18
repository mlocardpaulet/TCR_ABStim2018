psites <- c("Inpp5d_S1090", "Mapk3_T203+Y205", "Arhgef2_S955", "Zap70_Y491+Y492", "Zap70_T291", "Zap70_Y314", "Lat_T40", "Pstpip_S312", "Ubash3a_Y9")
source("RFunctions/PlotSite.R")
PlotSite(tab$psiteID[tab$GeneID %in% psites], tab)
tab$GeneID[tab$Regulation == "FALSE" & tab$GeneID %in% psites]

tab <- tab[tab$GeneID %in% psites,]

tab2 <- tab[, grepl("MeanLoops", names(tab))]
row.names(tab2) <- tab$psiteID


for (i in seq_len(nrow(tab))) {
  el <- tab$GeneID[i]
  sub <- tab2[i,]
  sub <- melt(sub)
  sub <- sub[!is.na(sub$value),]
  names(sub)[1] <-  "TimePoint"
  sub$TimePoint <- as.character(sub$TimePoint)
  sub$TimePoint <- gsub("MeanLoops_R[12345]_", "", sub$TimePoint)
  # I keep only the sites with a  minimum of 3 time points with points from a minimum of 2 biological repeats:
  #############################################
  temp <- table(sub$TimePoint)
  temp <- table(temp)
  if (sum(temp[names(temp) >= 2]) >= 3) {
    anovA <- aov(sub$value~sub$TimePoint, data = sub)
    anovB <- anova(anovA)
    pAnova <- anovB$"Pr(>F)"[1]
    Tuk <- TukeyHSD(anovA)
    pTuk <- extract_p(Tuk)[order(names(extract_p(Tuk)))]
    print(el)
    print(paste0("pAnova=", pAnova, " pTuk=", pTuk))
  } 
}

# With replacement of all the values:
# I take the 5% quantile of tab2 and a sd of 0.1
noise <- quantile(as.matrix(tab2), 0.05, na.rm = T)
noise <- rnorm(length(tab2[is.na(tab2)]), mean = noise, sd = 0.1)
hist(noise)

tab3 <- tab2
tab3[is.na(tab3)] <- noise
for (i in seq_len(nrow(tab))) {
  el <- tab$GeneID[i]
  if (grepl("_Y", el, fixed = T) | grepl("+Y", el, fixed = T)) {
    col <- grepl("pTyr", names(tab3))
  } else {
    col <- grepl("TiO2", names(tab3))
  }
  sub <- tab3[i,col]
  sub <- melt(sub)
  sub <- sub[!is.na(sub$value),]
  names(sub)[1] <-  "TimePoint"
  sub$TimePoint <- as.character(sub$TimePoint)
  sub$TimePoint <- gsub("MeanLoops_R[12345]_", "", sub$TimePoint)
  # I keep only the sites with a  minimum of 3 time points with points from a minimum of 2 biological repeats:
  #############################################
  temp <- table(sub$TimePoint)
  temp <- table(temp)
  if (sum(temp[names(temp) >= 2]) >= 3) {
    anovA <- aov(sub$value~sub$TimePoint, data = sub)
    anovB <- anova(anovA)
    pAnova <- anovB$"Pr(>F)"[1]
    Tuk <- TukeyHSD(anovA)
    pTuk <- extract_p(Tuk)[order(names(extract_p(Tuk)))]
    print(el)
    print(paste0("pAnova=", pAnova, " pTuk=", pTuk))
  } 
}

tab3$psiteID <- row.names(tab3)

PlotSite(tab$psiteID[tab$GeneID %in% psites[grepl("_Y", psites, fixed = T) | grepl("+Y", psites, fixed = T)]], tab3[,grepl("pTyr", names(tab3)) | names(tab3) == "psiteID"])
PlotSite(tab$psiteID[tab$GeneID %in% psites[!(grepl("_Y", psites, fixed = T) | grepl("+Y", psites, fixed = T))]], tab3[,grepl("TiO2", names(tab3)) | names(tab3) == "psiteID"])
