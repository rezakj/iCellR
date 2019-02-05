#' Add CITE-seq antibody-derived tags (ADT)
#'
#' This function takes a data frame of ADT values per cell and adds it to the iCellR object.
#' @param x An object of class iCellR.
#' @param adt.data A data frame containing ADT counts for cells.
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' my.obj <- add.adt(my.obj, adt.data = adt.data)
#' }
#'
#' @export
add.vdj <- function (x = NULL, vdj.data = "data.frame") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  if (class(vdj.data) != "data.frame") {
    stop("VDJ data should be a data frame object")
  }
  attributes(x)$vdj.data <- vdj.data
  return(x)
}
