#' Remove the cells that are in a cluster
#'
#' This function removes the cells from a designated cluster. Notice the cells will be removed from the main data (raw data would still have the original data).
#' @param x A data frame containing gene counts for cells.
#' @param clust.to.rm The name of the cluster to be removed.
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' my.obj <- clust.rm(my.obj, clust.to.rm = 5)
#' }
#'
#' @export
clust.rm <- function (x = NULL, clust.to.rm = "numeric") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  if (clust.to.rm == "numeric") {
    stop("you should choose a cluster number to remove")
  }
  # find cluster ids to remove
  DATA <- x@best.clust
  clustersToGo <- row.names(subset(DATA, DATA$clusters == clust.to.rm))
  # remove
  DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
  attributes(x)$best.clust <- DATA
  # main data
  DATA <- x@main.data
  DATA <- DATA[ , -which(names(DATA) %in% clustersToGo)]
  attributes(x)$main.data <- DATA
  # PCA
  DATA <- x@pca.data
  DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
  attributes(x)$pca.data <- DATA
  # tSNE
  DATA <- x@tsne.data
  DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
  attributes(x)$tsne.data <- DATA
  # umap
#  DATA <- x@umap.data
#  DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
#  attributes(x)$umap.data <- DATA
  #
  DATA <- x@tsne.data.3d
  DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
  attributes(x)$tsne.data.3d <- DATA
  return(x)
}
