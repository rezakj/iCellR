#' Pseudotime
#'
#' This function takes an object of class iCellR and marker genes for clusters and performs pseudotime analysis.
#' @param x An object of class iCellR.
#' @param marker.genes A list of marker genes for clusters.
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- pseudotime(my.obj, marker.genes = MyGenes, dim = 1:10)
#' }
#' @import gridExtra
#' @export
pseudotime <- function (x = NULL, marker.genes = "NULL", dim = 1:10) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  # geth the genes and scale them based on model
  if (dim(x@imputed.data)[1] == 0) {
    stop("data is not imputed")
  }
  DATA <- x@imputed.data
#  DATA <- x@main.data
  DATA <- subset(DATA, rownames(DATA) %in% marker.genes)
  data <- as.data.frame(scale(DATA))
  counts.pca <- prcomp(data, center = T, scale. = T)
  dataPCA = data.frame(counts.pca$rotation)
  DATA <- dataPCA
  DATA <- DATA[dim]
  myUMAP = umap(DATA, method = "umap-learn")
  myUMAP = as.data.frame((myUMAP$layout))
### return
attributes(x)$pca.data <- myUMAP
#attributes(x)$pca.data <- dataPCA
#cluster.plot(my.obj,plot.type = "pca",cell.color = "black",cell.transparency = 1,interactive = F)
#  attributes(x)$pseudo.mapA <- dataPCA
#  attributes(x)$pseudo.mapB <- myUMAP
  return(x)
}
