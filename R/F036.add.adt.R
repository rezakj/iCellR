#' Add CITE-seq antibody-derived tags (ADT)
#'
#' This function takes a data frame of ADT values per cell and adds it to the iCellR object.
#' @param x An object of class iCellR.
#' @param adt.data A data frame containing ADT counts for cells.
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' # Read your ADT data (in this case in ADT.tsv file) and add to your object as below
#' my.adt.data <- read.table("ADT.tsv")
#'
#' my.obj <- add.adt(my.obj, adt.data = my.adt.data)
#'
#' head(my.obj@adt.raw)[1:5]
#'
#' }
#'
#' @export
add.adt <- function (x = NULL, adt.data = "data.frame") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  if (class(adt.data) != "data.frame") {
    stop("ADT data should be a data frame object")
  }
  # Add ADT to gene names
  row.names(adt.data) <- paste("ADT", row.names(adt.data), sep = "_")
  row.names(adt.data) <- gsub("-",".",row.names(adt.data))
  row.names(adt.data) <- gsub("/","_or_",row.names(adt.data))
  row.names(adt.data) <- gsub(" ","_",row.names(adt.data))
  colnames(adt.data) <- gsub("-",".",colnames(adt.data))
  #
  adt.data[is.na(adt.data)] <- 0
  attributes(x)$adt.raw <- adt.data
  return(x)
}
