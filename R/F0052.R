#' Make BED Files
#'
#' This function takes peak marker files and makes the bed files per cluster.
#' @param x Peak marker file.
#' @return Bed files
#' @export
make.bed <- function (x = NULL) {
# get filed
  results <- subset(x,select=c(gene,clusters,AvExpInOtherClusters))
  peaks <- as.character(results$gene)
  peaks <- (gsub("\\.","_",peaks))
  peaks <- data.frame(do.call('rbind', strsplit(as.character(peaks),'_',fixed=TRUE)))
  results <- cbind(peaks,results)
#######
  My.clusters <- unique(results$clusters)
#######
  for(i in My.clusters){
    dat <- subset(results, results$clusters == i)
    dat <- subset(dat,select=c(X1,X2,X3,gene,AvExpInOtherClusters))
    MyName <- paste("peaks_cluster_",i,".bed",sep="")
    if(dim(dat)[1] > 0) {
      write.table(dat, MyName, sep="\t", quote = FALSE,row.names =FALSE, col.names = FALSE,)
    }
  }
}

