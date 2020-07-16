#' Choose top marker genes
#'
#' This function takes the marker genes info if chooses marker gene names for plots.
#' @param x An object of class iCellR.
#' @param topde Number of top differentially expressed genes to be choosen from each cluster, default = 10.
#' @param min.base.mean Minimum base mean of the genes to be chosen, default = 0.5.
#' @param filt.ambig Filter markers that are seen for more than one cluster, default = TRUE.
#' @param cluster Choose a cluster to find markers for. If 0, it would find markers for all clusters, , default = 0.
#' @return A set of gene names
#' @examples
#' marker.genes <- findMarkers(demo.obj,fold.change = 2,padjval = 0.1,uniq = TRUE)
#' top.markers(marker.genes, topde = 10, min.base.mean = 0.8)
#' @import Matrix
#' @export
top.markers <- function (x = NULL, topde = 10,
                         min.base.mean = 0.2,
                         filt.ambig = TRUE,
                         cluster = 0) {
  if (cluster == 0) {
    # get all clusters
  MyClusts <- (unique(x$clusters))
  # get rid of ambig genes (more than 1 cluster)
  data <- (as.data.frame(table(x$gene)))
  datanew <- (data[order(data$Freq, decreasing = TRUE),])
  datanew1 <- subset(datanew, datanew$Freq == 1)
  datanew1 <- as.character(datanew1$Var1)
  myDATA = x
  myDATA <- subset(myDATA, myDATA$gene %in% datanew1)
  if(filt.ambig == TRUE) {
    x = myDATA
  }
#  x <- x[order(x$baseMean,decreasing = T),]
  for (i in MyClusts) {
    DATA <- subset(x, x$clusters == i)
    DATA <- subset(DATA, DATA$baseMean >= min.base.mean)
    DATA <- as.character(head(DATA,topde)$gene)
    DatNmaes <- paste("mytopgenes", i, sep = "_")
    eval(call("<-", as.name(DatNmaes), DATA))
  }
  # cat them
  filenames <- ls(pattern="mytopgenes_")
  filenames <- filenames[order(nchar(filenames), filenames)]
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
