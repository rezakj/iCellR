#' Run UMAP on PCA data (Computes a manifold approximation and projection)
#'
#' This function takes an object of class iCellR and runs UMAP on PCA data.
#' @param x An object of class iCellR.
#' @param dims PC dimentions to be used for UMAP analysis.
#' @param method Character, implementation. Available methods are 'naive' (an implementation written in pure R) and 'umap-learn' (requires python package 'umap-learn'). Choose from "naive" and "umap-learn", defilt = "naive".
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- run.umap(my.obj, dims = 1:10)
#' }
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
  #  2 dimention
  #  if (clust.dim == 2) {
  # TransPosed <- t(TopNormLogScale)
  # DD <- phate(DATA)
  # phate.data <- as.data.frame(DD$embedding)
  #### method 2
  #devtools::install_github("jlmelville/uwot")
  # library(uwot)
  # myUMAP <- umap(DATA)
  # row.names(myUMAP) <- row.names(DATA)
  # myUMAP <- as.data.frame(myUMAP)
#  library(umap)
  myUMAP = umap(DATA, method = method)
  myUMAP = as.data.frame((myUMAP$layout))
  attributes(x)$umap.data <- myUMAP
# return
  return(x)
}
