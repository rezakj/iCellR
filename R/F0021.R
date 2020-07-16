#' Create a data frame of mean expression of genes per cluster
#'
#' This function takes an object of class iCellR and creates an average gene expression for every cluster.
#' @param x An object of class iCellR.
#' @param data.type Choose from "main" and "imputed", default = "main"
#' @param conds.to.avg Choose the conditions you want to average, default = NULL (all conditions).
#' @return An object of class iCellR.
#' @import progress
#' @examples
#' demo.obj <- clust.avg.exp(demo.obj)
#'
#' head(demo.obj@clust.avg)
#' @export
clust.avg.exp <- function (x = NULL,
                           data.type = "main",
                           conds.to.avg = NULL) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
      DATA <- x@best.clust
      MY.conds <- as.data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
      MY.conds <- as.character(as.matrix(MY.conds))
      DATA <- as.data.frame(cbind(DATA,MY.conds))
#######
      if (!is.null(conds.to.avg)) {
        DATA <- subset(DATA, DATA$MY.conds %in% conds.to.avg)
        message(paste(" "))
        message(paste("#############"))
        message(paste(" Averaging gene expression for all clusters in:",paste(conds.to.avg,collapse = ","),"..."))
        message(paste("#############"))
        message(paste(" "))
      }
###########
  # get data
  sampleCondition <- DATA$clusters
  conditions <- sort(unique(sampleCondition))
  DATA1 <- DATA
  ## get main data
  if (data.type == "main") {
    Table <- x@main.data
  }
  if (data.type == "imputed") {
    Table <- x@imputed.data
  }
#  Table = x@main.data
  for(i in conditions){
    IDs <- rownames(subset(DATA1, sampleCondition == i))
    DATA <- Table[ , which(names(Table) %in% IDs)]
    DATA <- as.matrix(DATA)
    message(paste(" Averaging gene expression for cluster:",i,"..."))
    message(paste("      Averaging",dim(DATA)[2],"cells ..."))
#    DATA <- Table[,row.names(subset(DATA1, sampleCondition == i))]
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
   message("All done!")
  return(x)
}
