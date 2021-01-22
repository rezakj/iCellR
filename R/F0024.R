#' Assign cluster number to cell ids
#'
#' This function takes an object of class iCellR and assigns cluster number to a vector of cell ids.
#' @param x An object of class iCellR.
#' @param my.gate A vector of cell ids.
#' @param to.clust A cluster id to be assigned to the provided cell ids.
#' @return An object of class iCellR.
#' @export
gate.to.clust <- function(x = NULL,
                        my.gate = NULL,
                        to.clust = 0) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  # geth data
  DATA <- (x@best.clust)
  gates <- readLines(my.gate)
#  DATA <- subset(DATA, row.names(DATA) %in% gates)
  DATA[which(row.names(DATA) %in% gates),1] <- to.clust
  #
  attributes(x)$best.clust <- DATA
  return(x)
}
