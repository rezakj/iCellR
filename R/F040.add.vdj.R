#' Add V(D)J recombination data
#'
#' This function takes a data frame of VDJ information per cell and adds it to the iCellR object.
#' @param x An object of class iCellR.
#' @param vdj.data A data frame containing VDJ information for cells.
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' Read your VDJ data (in this case in VDJ.tsv file) and add to your object as below
#'
#' my.vdj.data <- read.table("VDJ.tsv")
#'
#' VDJ <- prep.vdj(my.obj, adt.data = my.vdj.data)
#'
#' head(VDJ)
#'
#' my.obj <- add.vdj(my.obj, vdj.data = VDJ)
#'
#' head(my.obj@vdj.data)
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
