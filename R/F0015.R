#' Clustering the data
#'
#' This function takes an object of class iCellR and finds optimal number of clusters and clusters the data.
#' @param x An object of class iCellR.
#' @param clust.method the cluster analysis method to be used. This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans".
#' @param dist.method the distance measure to be used to compute the dissimilarity matrix. This must be one of: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski" or "NULL". By default, distance="euclidean". If the distance is "NULL", the dissimilarity matrix (diss) should be given by the user. If distance is not "NULL", the dissimilarity matrix should be "NULL".
#' @param index.method the index to be calculated. This should be one of : "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw", "all" (all indices except GAP, Gamma, Gplus and Tau), "alllong" (all indices with Gap, Gamma, Gplus and Tau included).
#' @param max.clust maximal number of clusters, between 2 and (number of objects - 1), greater or equal to min.nc.
#' @param min.clust minimum number of clusters, default = 2.
#' @param dims PCA dimentions to be use for clustering, default = 1:10.
#' @return An object of class iCellR.
#' @examples
#' demo.obj <- run.clustering(demo.obj,
#'                           clust.method = "kmeans",
#'                           dist.method = "euclidean",
#'                           index.method = "silhouette",
#'                           max.clust = 2,
#'                           min.clust = 2,
#'                           dims = 1:10)
#'
#'  head(demo.obj@best.clust)
#'
#' @import NbClust
#' @export
run.clustering <- function (x = NULL,
                          clust.method = "kmeans",
                          dist.method = "euclidean",
                          index.method = "silhouette",
                          max.clust = 25,
                          min.clust = 2,
                          dims = 1:10) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  #  cluster
  DATA <- (x@pca.data)[dims]
#  DATA <- x@umap.data
#  DATA <- x@tsne.data.3d
      nb <- NbClust(DATA,
                  distance = dist.method,
                  min.nc = min.clust,
                  max.nc = max.clust,
                  method = clust.method,
                  index = index.method)
    cluster.data <- nb
    attributes(x)$cluster.data <- cluster.data
    ####
    clustInfo = as.data.frame(x@cluster.data$Best.partition)
    colnames(clustInfo) <- "clusters"
    attributes(x)$best.clust <- clustInfo
  return(x)
}
