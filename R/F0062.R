#' Pseudotime
#'
#' This function takes an object of class iCellR and marker genes for clusters and performs pseudotime analysis.
#' @param x An object of class iCellR.
#' @param marker.genes A list of marker genes for clusters.
#' @param dims PC dimentions to be used, , default = 1:10.
#' @return An object of class iCellR.
#' @import gridExtra
#' @export
pseudotime <- function (x = NULL, marker.genes = "NULL", dims = 1:10) {
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
  counts.pca <- prcomp(data, center = TRUE, scale. = TRUE)
  dataPCA = data.frame(counts.pca$rotation)
  DATA <- dataPCA
  DATA <- DATA[dims]
  myUMAP = umap(DATA)
  myUMAP = as.data.frame((myUMAP))
### return
#attributes(x)$pca.data <- myUMAP
#attributes(x)$pca.data <- dataPCA
#cluster.plot(my.obj,plot.type = "pca",cell.color = "black",cell.transparency = 1,interactive = F)
  attributes(x)$pseudo.mapA <- dataPCA
  attributes(x)$pseudo.mapB <- myUMAP
  return(x)
}
