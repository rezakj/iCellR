#' Normalize ADT data.
#' This function takes data frame and Normalizes ADT data.
#' @param x An object of class iCellR.
#' @return An object of class iCellR
#' @export
norm.adt <- function (x = NULL) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  DATA <- x@adt.raw
  DATA <- as.data.frame(t(DATA))
  libSiz <- colSums(DATA)
  norm.facts <- as.numeric(libSiz) / mean(as.numeric(libSiz))
  dataMat <- as.matrix(DATA)
  normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  normalized <- as.data.frame(t(normalized))
#  DATA <- log2(DATA + 1)
  normalized[is.na(normalized)] <- 0
  DATA <- normalized
  attributes(x)$adt.main <- DATA
  return(x)
}
