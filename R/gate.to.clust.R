#' Run PCA on the main data
#'
#' This function takes an object of class iCellR and runs PCA on the main data.
#' @param x An object of class iCellR.
#' @param clust.method Choose from "base.mean.rank" or "gene.model", defult is "base.mean.rank".
#' @param top.rank A number taking the top genes ranked by base mean, defult = 500.
#' @param gene.list A list of genes to be used for PCA. If "clust.method" is set to "gene.model", defult = "my_model_genes.txt".
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- run.pca(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")
#' }
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
