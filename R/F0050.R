#' Sort and relabel the clusters randomly or based on pseudotime
#'
#' This function takes an object of class iCellR and re-ordersthe clusters based on pseudotime (distance).
#' @param x An object of class iCellR.
#' @param top.rank A number. Taking the top genes ranked by base mean, default = 500.
#' @param clust.method Choose from "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid", default = "complete".
#' @param dist.method Choose from "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski", default = "euclidean".
#' @param how.to.order Choose from "distance" and "random".
#' @return An object of class iCellR.
#' @import ape
#' @export
clust.ord <- function (x = NULL,
                       top.rank = 500,
                       dist.method = "euclidean",
                       clust.method = "complete",
                       how.to.order = "distance"
                       ) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  DATA <- x@best.clust
  if(!is.numeric(DATA$clusters)){
    stop("Clusters have to be numeric")
  }
  #####
  MyClust <- sort(unique(DATA$clusters))
  datC <- DATA$clusters + 1000
  for (i in 1:length(MyClust)) {
    From <- as.numeric(MyClust[i] + 1000)
    datC <- as.numeric(gsub(From,i,datC))
  }
  DATA$clusters <- datC
  # get data
  if(how.to.order == "distance") {
  x <- clust.avg.exp(x)
  data <- x@clust.avg
  row.names(data) <- data$gene
  data <- data[,-1]
  data <- data[ order(rowMeans(data), decreasing = TRUE), ]
  data <- head(data,top.rank)
  data <- dist(scale(t(data)), method = dist.method)
  hc <- hclust(data, method = clust.method)
  BestOrd <- hc$order
  ######
  orig <- DATA$clusters + 1000
  BestOrd.length <- c(1:length(BestOrd))
  for (i in BestOrd.length) {
    newClust <- BestOrd[i]
    newi <- newClust + 1000
    orig <- gsub(newi,i,orig)
  }
  DATA$clusters <- as.numeric(orig)
  x@best.clust <- DATA
  return(x)
  }
  if(how.to.order == "random") {
    MyClust <- sample(sort(unique(DATA$clusters)))
    datC <- DATA$clusters + 1000
    for (i in 1:length(MyClust)) {
      From <- as.numeric(MyClust[i] + 1000)
      datC <- as.numeric(gsub(From,i,datC))
    }
    DATA$clusters <- datC
    x@best.clust <- DATA
    return(x)
  }
}
