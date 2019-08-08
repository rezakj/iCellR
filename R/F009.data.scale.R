#' Scale data
#'
#' This function takes an object of class iCellR and scales the normalized data.
#' @param x An object of class iCellR.
#' @return An object of class iCellR.
#' @examples
#' my.obj <- data.scale(demo.obj)
#'
#' head(my.obj@scaled.data)[1:5]
#'
#' @export
data.scale <- function (x = NULL) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  DATA <- x@main.data
  NormLog = log2(DATA + 1)
  NormLog = as.data.frame(t(NormLog))
  TopNormLogScale <- as.data.frame(t(scale(NormLog, center=F)))
  attributes(x)$scaled.data <- TopNormLogScale
  return(x)
}

