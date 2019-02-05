#' Run UMAP on PCA data
#'
#' This function takes an object of class iCellR and runs UMAP on PCA data.
#' @param x An object of class iCellR.
#' @param dims PC dimentions to be used for UMAP analysis.
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- run.magic(my.obj)
#' }
#' @import umap
#' @export
run.magic <- function (x = NULL) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
### load packages
  # https://github.com/KrishnaswamyLab/MAGIC/tree/master/Rmagic
  # https://www.analyticsvidhya.com/blog/2016/03/tutorial-powerful-packages-imputing-missing-values/
  require(Rmagic)
  require(viridis)
  require(phateR)
  # get data
  DATA <- x@main.data
  DATA <- as.data.frame(t(DATA))
  data_MAGIC <- magic(DATA, genes="all_genes")
  DATA <- as.data.frame(t(data_MAGIC$result))
  # return
  attributes(x)$main.data <- DATA
  # return
  return(x)
}
