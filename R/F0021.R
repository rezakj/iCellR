#' Create a data frame of mean expression of genes per cluster
#'
#' This function takes an object of class iCellR and creates an average gene expression for every cluster.
#' @param x An object of class iCellR.
#' @param data.type Choose from "main", "atac", "atac.imputed" and "imputed", default = "main"
#' @param conds.to.avg Choose the conditions you want to average, default = NULL (all conditions).
#' @param rounding.digits integer indicating the number of decimal places (round) or significant digits (signif) to be used.
#' @param round.num Rounding of Numbers, default = FALSE.
#' @return An object of class iCellR.
#' @import progress
#' @export
clust.avg.exp <- function (x = NULL,
                           data.type = "main",
                           conds.to.avg = NULL,
                           rounding.digits = 4,
                           round.num = FALSE) {
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
  if (data.type == "atac") {
    Table <- x@atac.main
  }
  if (data.type == "atac.imputed") {
    Table <- x@atac.imputed
  }
#  Table = x@main.data
  datalist <- list()
  ###
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
#    datalist[[i]] <- DATA
  }
#  multmerge = function(mypath){
#    filenames=list.files(pattern="meanExp")
#    datalist = lapply(filenames, function(x){read.table(file=x,header=T)})
#    Reduce(function(x,y) {merge(x,y)}, datalist)
#  }
   filenames <- ls(pattern="cluster_")
   filenames <- filenames[order(nchar(filenames))]
   datalist <- mget(filenames)
   MeanExpForClusters <- Reduce(function(x,y) {merge(x,y)}, datalist)
#
#  MeanExpForClusters <- multmerge()
#  file.remove(list.files(pattern="meanExp"))
#   MyMat <- as.matrix(MeanExpForClusters)
#   NameS <- colnames(MeanExpForClusters)
   #######
#   MeanExpForClusters <- MyMat[,order(nchar(NameS),NameS)]
#   MeanExpForClusters <-as.data.frame(MeanExpForClusters)
   ##############
   data <- MeanExpForClusters
   row.names(data) <- data$gene
   data <- data[,-1]
#   data = as.data.frame(sapply(data, as.numeric))
   if (round.num == TRUE) {
     data <- round(data, digits = rounding.digits)
   }
   data <- cbind(gene=MeanExpForClusters$gene,data)
#   MeanExpForClusters <- MeanExpForClusters[order(nchar(colnames(MeanExpForClusters)),colnames(MeanExpForClusters))]
      attributes(x)$clust.avg <- data
   message("All done!")
  return(x)
}
