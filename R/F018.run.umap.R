#' Run UMAP on PCA Data (Computes a manifold approximation and projection)
#'
#' This function takes an object of class iCellR and runs UMAP on PCA data.
#' @param x An object of class iCellR.
#' @param dims PC dimentions to be used for UMAP analysis.
#' @param method Character, implementation. Available methods are 'naive' (an implementation written in pure R) and 'umap-learn' (requires python package 'umap-learn'). Choose from "naive" and "umap-learn", default = "naive".
#' @return An object of class iCellR.
#' @examples
#' demo.obj <- run.umap(demo.obj, dims = 1:10)
#' head(demo.obj@umap.data)
#'
#' @import umap
#' @export
run.umap <- function (x = NULL,
                         dims = 1:10,
                      method = "naive") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  # https://github.com/lmcinnes/umap
  # get PCA data
  DATA <- x@pca.data
  DATA <- DATA[dims]
#  data <- as.matrix(DATA)
#  Rphenograph_out <- Rphenograph(data, k = 45)
#  modularity(Rphenograph_out[[2]])
#  MyClusts = membership(Rphenograph_out[[2]])
#  MyClusts <- as.data.frame(as.matrix(MyClusts))
#  row.names(MyClusts) <- row.names(data)
#  colnames(MyClusts) <- "clusters"
#  my.obj@best.clust <- MyClusts
#  head(my.obj@best.clust)
  ##########
#  DATA <- my.obj@imputed.data
#  DATA <- my.obj@main.data
#  raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = TRUE), ]
#  topGenes <- head(raw.data.order,500)
#  DATA<- as.data.frame(t(as.data.frame(scale(topGenes))))
#  head(DATA)[1:5]
#  myUMAP = umap(DATA, method = "umap-learn")
#  myUMAP = as.data.frame((myUMAP$layout))
#  my.obj@umap.data <- myUMAP
#  head(my.obj@umap.data)
#  library(umap)
  myUMAP = umap(DATA, method = method)
  myUMAP = as.data.frame((myUMAP$layout))
  attributes(x)$umap.data <- myUMAP
# return
  return(x)
}
