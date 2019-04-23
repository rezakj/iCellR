#' Change the cluster number or re-name them
#'
#' This function re-names the clusters in the best.clust slot of the iCellR object.
#' @param x An object of class iCellR.
#' @param change.clust The name of the cluster to be changed.
#' @param to.clust The new name for the cluster.
#' @param clust.reset Reset to the original clustering.
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- change.clust(my.obj, change.clust = 3, to.clust = 1)
#' my.obj <- change.clust(my.obj, change.clust = 2, to.clust = "B Cell")
#' my.obj <- change.clust(my.obj, clust.reset = T)
#' }
#' @export
change.clust <- function (x = NULL,
                          change.clust = 0,
                          to.clust = 0,
                          clust.reset = F) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
###################
  if (clust.reset == F) {
    DATA <- (x@best.clust)
    DATA[DATA == change.clust] <- to.clust
    attributes(x)$best.clust <- DATA
  }
# reset
  if (clust.reset == T) {
    clustInfo = as.data.frame(x@cluster.data$Best.partition)
    colnames(clustInfo) <- "clusters"
    attributes(x)$best.clust <- clustInfo
  }
##############
  return(x)
}
