#' Merge RNA and ADT data
#'
#' This function is to merge the RNA and ADT data to the main.data slot of the iCellR object.
#' @param x An object of class iCellR.
#' @param adt.data Choose from raw or main (normalized) ADT data, default = "raw".
#' @return An object of class iCellR
#' @export
adt.rna.merge <- function (x = NULL, adt.data = "raw") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ### get ADT data
  if (adt.data == "raw") {
  DATA <- as.data.frame(t(x@adt.raw))
  }
  if (adt.data == "main") {
  DATA <- as.data.frame(t(x@adt.main))
  }
  ### get RNA data
  MainData <- x@main.data
  RmADTs <- grep("^ADT_", row.names(MainData), value = TRUE)
  MainData <- subset(MainData,! row.names(MainData) %in% RmADTs)
  DATArna <- as.data.frame(t(MainData))
  ### merge
  merged.data <- merge(DATA, DATArna, by="row.names", all.x=FALSE, all.y=TRUE)
  merged.data[is.na(merged.data)] <- 0
  rownames(merged.data) <- merged.data$Row.names
  merged.data <- merged.data[,-1]
  merged.data <- as.data.frame(t(merged.data))
  attributes(x)$main.data <- merged.data
  return(x)
}
