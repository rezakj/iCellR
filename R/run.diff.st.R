#' Differentiation spacetime
#'
#' This function takes an object of class iCellR and finds their progenitor cells based on distance.
#' @param x An object of class iCellR.
#' @param dist.method the distance measure to be used to compute the dissimilarity matrix. This must be one of: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski" or "NULL". By default, distance="euclidean". If the distance is "NULL", the dissimilarity matrix (diss) should be given by the user. If distance is not "NULL", the dissimilarity matrix should be "NULL".
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- dst(my.obj)
#' }
#' @import NbClust
#' @export
run.diff.st <- function (x = NULL,
                         dist.method = "euclidean") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  #  cluster
  DATA <- x@tsne.data
  data <- data.matrix(DATA)
  data <- dist(data, method = dist.method)
  data <- as.data.frame(as.matrix(data))[1]
#  data <-  as.data.frame(log2(as.matrix(data) + 1))
  DATA <- x@tsne.data
  data <- cbind(data, DATA)
#  My.Clust <- x@best.clust
 # data <- cbind(data, My.Clust)
  colnames(data) <- c("differentiation.time", "spatial.dim.1","spatial.dim.2")
  # plot
#  col.legend <- as.factor(data$clusters)
#ggplot(data,aes(x=differentiation.time,y=spatial.dim.1, col=spatial.dim.2)) +
#  geom_point()
#ggplot(data,aes(x=differentiation.time,y=spatial.dim.1, col=col.legend)) +
#  geom_point()
# + scale_y_continuous(trans = "log1p") +
#  scale_x_continuous(trans = "log1p")
#DATA = data
#plot_ly(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3],
#        color = col.legend, opacity = cell.transparency,
#        marker = list(size = cell.size + 2))
  # return
  attributes(x)$diff.st.data <- data
    return(x)
}
