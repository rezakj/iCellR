#' Run tSNE on PCA data
#'
#' This function takes an object of class iCellR and runs tSNE on PCA data.
#' @param x An object of class iCellR.
#' @param dims PC dimentions to be used for tSNE analysis.
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- run.pc.tsne(my.obj, dims = 1:10)
#' }
#' @import Rtsne
#' @export
run.pc.tsne <- function (x = NULL,
                      dims = 1:10) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
# get PCA data
  DATA <- x@pca.data
  TransPosed <- DATA[dims]
  #  2 dimention
  #  if (clust.dim == 2) {
  # TransPosed <- t(TopNormLogScale)
  tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = 2)
  tsne.data = as.data.frame(tsne$Y)
  tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
  rownames(tsne.data) <- tsne.data$cells
  tsne.data <- tsne.data[,-1]
  attributes(x)$tsne.data <- tsne.data
  #  }
  # choose 3 demention
  # tSNE
  # TransPosed <- t(TopNormLogScale)
  tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = 3)
  tsne.data = as.data.frame(tsne$Y)
  tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
  rownames(tsne.data) <- tsne.data$cells
  tsne.data <- tsne.data[,-1]
  attributes(x)$tsne.data.3d <- tsne.data
  return(x)
}
