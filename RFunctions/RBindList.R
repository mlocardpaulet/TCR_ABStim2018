RBindList <- function(list) {
  # combine tables of a list using rbind
  tab <- rbind(list[[1]], list[[2]])
  if (length(list)>2) {
    for (i in 3:length(list)) {
      tab <- rbind(tab, list[[i]])
    }
  }
  return(tab)
}