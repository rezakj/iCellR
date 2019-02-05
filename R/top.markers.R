#' Choose top marker genes
#'
#' This function takes the marker genes info if chooses marker gene names for plots.
#' @param x An object of class iCellR.
#' @param topde Number of top differentialy expressed genes to be choosen from each cluster, defult = 10.
#' @param min.base.mean Minimum base mean of the genes to be choosen, defult = 0.5.
#' @return A set of gene names
#' @examples
#' \dontrun{
#' MyGenes <- top.markers(marker.genes, topde = 10, min.base.mean = 0.8)
#' }
#' @import Matrix
#' @export
top.markers <- function (x = NULL, topde = 10, min.base.mean = 0.2, cluster = 0) {
  if (cluster == 0) {
  MyClusts <- (unique(x$clusters))
#  x <- x[order(x$baseMean,decreasing = T),]
  for (i in MyClusts) {
    DATA <- subset(x, x$clusters == i)
    DATA <- subset(DATA, DATA$baseMean >= min.base.mean)
    DATA <- as.character(head(DATA,topde)$gene)
    DatNmaes=paste("topgenes",i,sep="_")
    eval(call("<-", as.name(DatNmaes), DATA))
  }
  # cat them
  filenames <- ls(pattern="topgenes_")
  datalist <- mget(filenames)
  topGenes <- as.character(do.call("c", datalist))
  }
  if (cluster != 0) {
    MyClusts <- cluster
    DATA <- subset(x, x$clusters == MyClusts)
    DATA <- subset(DATA, DATA$baseMean >= min.base.mean)
    topGenes <- as.character(head(DATA,topde)$gene)
  }
    return(topGenes)
  }


