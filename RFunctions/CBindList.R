CBindList <- function(l) {
  # combine tables of a list l using cbind
  if (length(l)==1) {tab <- l[[1]]} else {
    if (class(l[[1]])!="matrix" | class(l[[1]])!="data.frame") {
      tab <- cbind(l[[1]], l[[2]])
      if (length(l)>2) {
        for (i in 3:length(l)) {
          tab <- cbind(tab, l[[i]])
        }
      } 
    } else {
      l <- l[2:length(l)]
      CBindList(l)
    }
  }
  return(tab)
}