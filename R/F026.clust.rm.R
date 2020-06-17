#' Remove the cells that are in a cluster
#'
#' This function removes the cells from a designated cluster. Notice the cells will be removed from the main data (raw data would still have the original data).
#' @param x A data frame containing gene counts for cells.
#' @param clust.to.rm The name of the cluster to be removed.
#' @return An object of class iCellR
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
  if (dim(x@main.data)[1] != 0) {
    DATA <- x@main.data
    DATA <- DATA[ , -which(names(DATA) %in% clustersToGo)]
    attributes(x)$main.data <- DATA
  }
  #
  # imputed data
  if (dim(x@imputed.data)[1] != 0) {
    DATA <- x@imputed.data
    DATA <- DATA[ , -which(names(DATA) %in% clustersToGo)]
    attributes(x)$imputed.data <- DATA
  }
  # PCA
  if (dim(x@pca.data)[1] != 0) {
    DATA <- x@pca.data
    DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
    attributes(x)$pca.data <- DATA
  }
  # tSNE
  if (dim(x@tsne.data)[1] != 0) {
    DATA <- x@tsne.data
    DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
    attributes(x)$tsne.data <- DATA
  }
  # umap
  if (dim(x@umap.data)[1] != 0) {
    DATA <- x@umap.data
    DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
    attributes(x)$umap.data <- DATA
  }
  #
  # knetl 2D
  if (dim(x@knetl.data)[1] != 0) {
    DATA <- x@knetl.data
    DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
    attributes(x)$knetl.data <- DATA
  }
  # knetl 3D
  if (dim(x@knetl.data.3d)[1] != 0) {
    DATA <- x@knetl.data.3d
    DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
    attributes(x)$knetl.data.3d <- DATA
  }
  # tsne 3d
  if (dim(x@tsne.data.3d)[1] != 0) {
    DATA <- x@tsne.data.3d
    DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
    attributes(x)$tsne.data.3d <- DATA
  }
  message("The cells in the cluster of your choice are removed from: main.data, imputed.data , pca.data, tsne.data, tsne.data.3d, knetl.data, knetl.data.3d and umap.data.")
  return(x)
}
