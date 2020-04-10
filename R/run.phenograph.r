#' Clustering the data
#'
#' This function takes an object of class iCellR and finds optimal number of clusters and clusters the data.
#' @param x An object of class iCellR.
#' @param k integer; number of nearest neighbours (default:45)
#' @param dims PCA dimentions to be use for clustering, default = 1:10.
#' @return An object of class iCellR.
#' @export
run.phenograph <- function (x = NULL,
                            k = 45,
                            dims = 1:10) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  #  cluster
  DATA <- (x@pca.data)[dims]
  #  DATA <- x@umap.data
  #  DATA <- x@tsne.data.3d
    ######
  ######
      data <- as.matrix(DATA)
      Rphenograph_out <- Rphenograph(data, k = k)
      MyClusts = membership(Rphenograph_out[[2]])
      MyClusts <- as.data.frame(as.matrix(MyClusts))
      row.names(MyClusts) <- row.names(data)
      colnames(MyClusts) <- "clusters"
  attributes(x)$best.clust <- MyClusts
  return(x)
}
