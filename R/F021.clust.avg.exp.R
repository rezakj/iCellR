#' Create a data frame of mean expression of genes per cluster
#'
#' This function takes an object of class iCellR and creates an average gene expression for every cluster.
#' @param x An object of class iCellR.
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- clust.avg.exp(my.obj)
#' }
#' @export
clust.avg.exp <- function (x = NULL) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
      DATA <- x@best.clust
  # get data
  sampleCondition <- DATA$clusters
  conditions <- unique(sampleCondition)
  DATA1 <- DATA
  Table = data.matrix(x@main.data)
  for(i in conditions){
    DATA <- Table[,row.names(subset(DATA1, sampleCondition == i))]
    DATA <- apply(DATA, 1, function(DATA) {mean(DATA)})
    DATA <- as.data.frame(DATA)
    Name=paste("meanExp_cluster",i,".txt",sep="_")
    NameCol=paste("cluster",i,sep="_")
    colnames(DATA) <- NameCol
    DATA <- cbind(gene = rownames(DATA), DATA)
    rownames(DATA) <- NULL
    eval(call("<-", as.name(NameCol), DATA))
#    head(DATA)
#    write.table((DATA),file=Name,sep="\t", row.names =F)
  }
#  multmerge = function(mypath){
#    filenames=list.files(pattern="meanExp")
#    datalist = lapply(filenames, function(x){read.table(file=x,header=T)})
#    Reduce(function(x,y) {merge(x,y)}, datalist)
#  }
   filenames <- ls(pattern="cluster_")
   datalist <- mget(filenames)
   MeanExpForClusters <- Reduce(function(x,y) {merge(x,y)}, datalist)
#
#  MeanExpForClusters <- multmerge()
#  file.remove(list.files(pattern="meanExp"))
   MeanExpForClusters <- MeanExpForClusters[order(nchar(colnames(MeanExpForClusters)),colnames(MeanExpForClusters))]
   attributes(x)$clust.avg <- MeanExpForClusters
  return(x)
}
