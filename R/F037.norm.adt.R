#' Normalize ADT data.
#' This function takes data frame and Normalizes ADT data.
#' @param x An object of class iCellR.
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' my.obj <- make.obj(my.obj)
#' }
#'
#' @export
norm.adt <- function (x = NULL) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  DATA <- x@adt.raw
  DATA <- log2(DATA + 1)
  attributes(x)$adt.main <- DATA
  return(x)
}
