#' Find optimal number of clusters.
#'
#' This function takes an object of class iCellR and finds optimal number of clusters based on three methods.
#' @param x An object of class iCellR.
#' @param max.clust Maximum number of clusters allowed to optimize, defult = 20.
#' @param gap.stat.nboot Number of boots, defult = 100 (500 is desired but slower).
#' @param verbose If TRUE will show the processing.
#' @param clust.type Choose from "tsne","pca" or "distance", defult = "tsne".
#' @param opt.method Choose from "elbow.wss", "silhouette", "gap.stat" or "", defult = "all".
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' opt.clust.num(my.obj, max.clust = 10, clust.type = "tsne", opt.method = "silhouette")
#' }
#' @import gridExtra
#' @export
opt.clust.num <- function (x = NULL,
                           max.clust = 20,
                           gap.stat.nboot = 100,
                           verbose = TRUE,
                           clust.type = "tsne",
                           opt.method = "silhouette") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
# get data
    if (clust.type == "tsne") {
      df <- (x@tsne.data)[1:2]
    }
    if (clust.type == "pca") {
      df <- (x@pca.data)[1:2]
    }
  #
  if (clust.type == "distance") {
    df <- x@dist.data
  }
  if (opt.method == "all") {
    Elbow = fviz_nbclust(df, kmeans, method = "wss", k.max = max.clust) +
      geom_vline(xintercept = 4, linetype = 2) +
      labs(subtitle = "Elbow method")
    Silhouette = fviz_nbclust(df, kmeans, method = "silhouette", k.max = max.clust) +
      labs(subtitle = "Silhouette method")
    set.seed(123)
    Gap.statistic = fviz_nbclust(df, kmeans,
                                 nstart = 25,
                                 method = "gap_stat",
                                 nboot = gap.stat.nboot,
                                 k.max = max.clust) +
      labs(subtitle = "Gap statistic method")
    return(grid.arrange(Elbow, Silhouette,Gap.statistic, nrow = 3))
  }
  if (opt.method == "elbow.wss") {
    Elbow = fviz_nbclust(df, kmeans, method = "wss", k.max = max.clust) +
      geom_vline(xintercept = 4, linetype = 2) +
      labs(subtitle = "Elbow method")
    return(Elbow)
  }
  if (opt.method == "silhouette") {
      geom_vline(xintercept = 4, linetype = 2)
      Silhouette = fviz_nbclust(df, kmeans, method = "silhouette", k.max = max.clust) +
      labs(subtitle = "Silhouette method")
      return(Silhouette)
  }
  if (opt.method == "gap.stat") {
    geom_vline(xintercept = 4, linetype = 2)
    Gap.statistic = fviz_nbclust(df, kmeans,
                                 nstart = 25,
                                 method = "gap_stat",
                                 nboot = gap.stat.nboot,
                                 k.max = max.clust) +
      labs(subtitle = "Gap statistic method")
    return(Gap.statistic)
  }
}
