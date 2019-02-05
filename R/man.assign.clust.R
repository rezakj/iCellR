#' Hierarchical clustering based on tSNE
#'
#' This function takes an object of class iCellR and performs hierarchical clustering based on tSNE data.
#' @param x An object of class iCellR.
#' @param clust.num Number of clusters to be made.
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' assign.clust(my.obj, clust.num = 7)
#' }
#' @export
man.assign.clust <- function (x = NULL,
                          clust.num = 0) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  if ( clust.num == 0) {
    stop("please provite the optimal number of clusters")
  }
  ClustNum = clust.num
#####
    DATA <- (x@tsne.data)
    fit_cluster_hierarchical = hclust(dist(scale(DATA)))
    DATA$clusters = factor(cutree(fit_cluster_hierarchical, k=ClustNum))
    DATAclusters = as.data.frame(cutree(fit_cluster_hierarchical, k = ClustNum))
    colnames(DATAclusters) <- c("clusters")
    attributes(x)$best.clust <- DATAclusters
  #
  return(x)
}
