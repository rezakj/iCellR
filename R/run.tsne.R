#' Run tSNE on the main data
#'
#' This function takes an object of class iCellR and runs tSNE on the main data.
#' @param x An object of class iCellR.
#' @param clust.method Choose from "base.mean.rank" or "gene.model", defult is "base.mean.rank".
#' @param top.rank A number taking the top genes ranked by base mean, defult = 500.
#' @param gene.list A list of genes to be used for tSNE analysis. If "clust.method" is set to "gene.model", defult = "my_model_genes.txt".
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- run.tsne(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")
#' }
#' @import Rtsne
#' @export
run.tsne <- function (x = NULL,
                      clust.method = "base.mean.rank",
                      top.rank = 500,
                      gene.list = "character") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
#  if (clust.dim != 2 && clust.dim != 3) {
#    stop("clust.dim should be either 2 or 3")
#  }
  if (clust.method == "dispersed.genes" && clust.method == "both") {
    stop("dispersed.genes and both are not implemented yet")
  }
  # geth the genes and scale them based on model
  DATA <- x@main.data
  # model base mean rank
  if (clust.method == "base.mean.rank") {
    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = T), ]
    topGenes <- head(raw.data.order,top.rank)
    TopNormLogScale <- log(topGenes + 0.1)
#    TopNormLogScale <- t(TopNormLogScale)
#    TopNormLogScale <- as.data.frame(t(scale(TopNormLogScale)))
  }
  # gene model
  if (clust.method == "gene.model") {
    if (gene.list == "character") {
      stop("please provide gene names for clustering")
    } else {
      genesForClustering <-gene.list
      topGenes <- subset(DATA, rownames(DATA) %in% genesForClustering)
      TopNormLogScale <- log(topGenes + 0.1)
#      TopNormLogScale <- t(TopNormLogScale)
#      TopNormLogScale <- as.data.frame(t(scale(TopNormLogScale)))
    }
  }
#  2 dimention
#  if (clust.dim == 2) {
      TransPosed <- t(TopNormLogScale)
      tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = 2)
      tsne.data = as.data.frame(tsne$Y)
      tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
      rownames(tsne.data) <- tsne.data$cells
      tsne.data <- tsne.data[,-1]
      attributes(x)$tsne.data <- tsne.data
#  }
# choose 3 demention
  # tSNE
      TransPosed <- t(TopNormLogScale)
      tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = 3)
      tsne.data = as.data.frame(tsne$Y)
      tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
      rownames(tsne.data) <- tsne.data$cells
      tsne.data <- tsne.data[,-1]
      attributes(x)$tsne.data.3d <- tsne.data
  return(x)
}
