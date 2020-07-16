#' Add V(D)J recombination data
#'
#' This function takes a data frame of VDJ information per cell and adds it to the iCellR object.
#' @param x An object of class iCellR.
#' @param vdj.data A data frame containing VDJ information for cells.
#' @return An object of class iCellR
#' @examples
#' my.vdj <- read.csv(file = system.file('extdata', 'all_contig_annotations.csv',
#'           package = 'iCellR'),
#'           as.is = TRUE)
#' head(my.vdj)
#' dim(my.vdj)
#'
#' My.VDJ <- prep.vdj(vdj.data = my.vdj, cond.name = "NULL")
#' head(My.VDJ)
#' dim(My.VDJ)
#'
#' my.obj <- add.vdj(demo.obj, vdj.data = My.VDJ)
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
