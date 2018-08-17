# Modified funtions without using the column Gene name but the Fasta headers instead:

MultiPhosphoParsingMQ <- function(tab) {
  # Input a Phospho(STY) table from MaxQuant and return the same table with additional rows containing info. for the multiply phosphorylated peptides (quan. values corresponding to several phosphorylation sites close together)
  
  # Remove the "_" in the Proteins field if present:
  tab$Proteins <- gsub("_", "|", tab$Proteins, fixed = T)
  Accession <- sapply(tab$Proteins, function(x) {
    strsplit(x, ";", fixed = T)
  })
  Accession <- lapply(Accession, function(x) {
    sapply(x, function(y) {
      strsplit(y, "|", fixed = T)[[1]][2]
    })
  })
  Accession <- lapply(Accession, function(x) {
    paste(x, collapse = ";")
  })
  tab$Accession <- unlist(Accession)
  
  # Add columns at the end of the table with the intensities derived from all the mono and multiply phosphorylated peptides:
  mat <- tab[,grepl("___1", names(tab))]
  names(mat) <- gsub("___1", "_Parsed", names(mat))
  nc <- ncol(mat) # Keep the number of columns to amend
  tab$phosphoSites <- paste0(tab$Accession, "_", tab$'Amino.acid', tab$'Positions.within.proteins')
  tab <- cbind(tab, mat) # Add the values of the mono-phosphorylated peptides.
  
  # Get the multiply phosphorylated values:
  mat <- tab[,c(which(grepl("___[23456789]", names(tab))), which(names(tab) == "Accession"), which(names(tab) == "Amino.acid"), which(names(tab) == "Positions.within.proteins"))]
  mat <- mat[mat$Accession != "", ]
  vecphosnum <- sapply(names(mat), function(x) {
    as.numeric(strsplit(x, "___", fixed = T)[[1]][2])
  })
  vecphosnum <- vecphosnum[!is.na(vecphosnum)]
  for (j in sort(unique(vecphosnum))) {
    matj <- mat[,c(which(vecphosnum == j), ncol(mat)-2, ncol(mat)-1, ncol(mat))]
    matj[,1:(ncol(matj)-3)] <- apply(matj[,1:(ncol(matj)-3)], 2, as.numeric)
    if (sum(matj[,1:(ncol(matj)-3)], na.rm = T) > 0) {
      matj <- matj[rowSums(matj[,1:(ncol(matj)-3)], na.rm = T) > 0,]
      matj <- matj[order(matj$'Positions.within.proteins'),]
      matj <- matj[order(matj$Accession),]
      matj[matj == 0] <- NA
      start <- 1
      for (el in unique(as.character(matj$Accession))) {
        matel <- matj[matj$Accession == el,]
        matel$sites <- paste0(matel$'Amino.acid', matel$'Positions.within.proteins')
        for (i in seq_len((ncol(matel)-4))) {
          veci <- matel[,i]
          veci <- veci[!is.na(veci)]
          tveci <- table(veci)
          tveci <- tveci[tveci == j]
          if (!is.null(tveci)) {
            for (k in seq_along(tveci)) {
              int <- names(tveci)[k]
              psiteID <- paste0(el, "_", paste(as.character(matel$sites[which(matel[,i] == int)]), collapse = "+"))
              if (start == 1) {
                mappend <- matrix(c(psiteID, i, int), nrow = 1)
                colnames(mappend) <- c("psiteID", "col", "intensity")
                start <- 0
              } else {
                mappend <- rbind(mappend, c(psiteID, i, int))
              }
            }
          }
        }
      }
      # Create matrix of new data:
      nr <- length(unique(mappend[,1]))
      mat2 <- matrix(nrow = nr, ncol = nc)
      row.names(mat2) <- unique(mappend[,1])
      for (k in seq_along(row.names(mat2))) {
        el <- row.names(mat2)[k]
        matel <- mappend[mappend[,1] == el,]
        if (class(matel) == "character") {
          mat2[k,as.numeric(matel[2])] <- as.numeric(matel[3])
        } else if (class(matel) == "matrix"){
          for (i in seq_len(nrow(matel))) {
            mat2[k,as.numeric(matel[i,2])] <- as.numeric(matel[i,3])
          }
        }
      }
      mat2 <- cbind(row.names(mat2), mat2)
      mat2 <- cbind(matrix(ncol = ncol(tab)-ncol(mat2), nrow = nrow(mat2)), mat2)
      colnames(mat2) <- colnames(tab)
      row.names(mat2) <- NULL
      tab <- rbind(tab, mat2)
    }
  }
  # Add more info the the added rows:
  tab$Accession <- as.character(tab$Accession)
  tab$Accession[is.na(tab$Accession)] <- sapply(as.character(tab$phosphoSites[is.na(tab$Accession)]) , function(x) {
    strsplit(x, "_", fixed = T)[[1]][1]
  })
  tab$'Gene.names' <- sapply(as.character(tab$`Fasta.headers`), function(x) {
    strsplit(x, "GN=", fixed = T)[[1]][2]
  })
  tab$'Gene.names' <- sapply(tab$`Gene.names`, function(x) {
    strsplit(x, " ", fixed = T)[[1]][1]
  })
  if (sum(names(tab) == "Gene.names") == 0 & sum(names(tab) == "Gene.Name") == 1) {
    names(tab)[names(tab) == "Gene.Name"] <- "Gene.names"
  }
  mappingGenes <- tab[,c("Accession", "Gene.names", "Fasta.headers")]
  mappingGenes <- mappingGenes[!duplicated(mappingGenes),]
  mappingGenes <- mappingGenes[!is.na(mappingGenes$Accession),]
  mappingGenes <- mappingGenes[!is.na(mappingGenes$Gene.names) & !is.na(mappingGenes$Gene.names),]
  tab$'Gene.names' <- as.character(tab$'Gene.names')
  tab$'Gene.names'[is.na(tab$'Gene.names')] <- as.character(mappingGenes$'Gene.names'[match(tab$Accession[is.na(tab$'Gene.names')], mappingGenes$Accession)])
  tab$'Fasta.headers' <- as.character(tab$'Fasta.headers')
  tab$'Fasta.headers'[is.na(tab$'Fasta.headers')] <- as.character(mappingGenes$'Fasta.headers'[match(tab$Accession[is.na(tab$'Fasta.headers')], mappingGenes$Accession)])
  return(tab)
}
