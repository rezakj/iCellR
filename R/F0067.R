#' Clustering the data
#'
#' This function takes an object of class iCellR and finds optimal number of clusters and clusters the data.
#' @param x An object of class iCellR.
#' @param k integer; number of nearest neighbours (default:45)
#' @param data.type Choose between "tsne", "pca", "umap", default = "pca".
#' @param dims PCA dimentions to be use for clustering, default = 1:10.
#' @return An object of class iCellR.
#' @export
run.phenograph <- function (x = NULL,
                            k = 100,
                            data.type = "pca",
                            dims = 1:10) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  #  cluster
  start_time1 <- Sys.time()
  #  cluster
  if(data.type == "pca") {
    DATA <- (x@pca.data)[dims]
    message(paste("Clustering based on PCA data"))
  }
  if(data.type == "umap") {
    DATA <- (x@umap.data[1:2])
    message(paste("Clustering based on UMAP data"))
  }
  if(data.type == "tsne") {
    DATA <- (x@tsne.data)
    message(paste("Clustering based on tSNE data"))
  }
  if(data.type == "knetl") {
    DATA <- (x@knetl.data)
    message(paste("Clustering based on KNetL data"))
  }
  ######
      data <- as.matrix(DATA)
      Rphenograph_out <- Rphenograph(data, k = k)
      MyClusts = membership(Rphenograph_out[[2]])
      MyClusts <- as.data.frame(as.matrix(MyClusts))
      row.names(MyClusts) <- row.names(data)
      colnames(MyClusts) <- "clusters"
  attributes(x)$best.clust <- MyClusts
  end_time1 <- Sys.time()
  Time = difftime(end_time1,start_time1,units = "mins")
  Time = round(as.numeric(Time),digits = 2)
  message(paste(" "))
  message(paste("Total time",Time,"mins"))
  message(paste("All done!"))
  return(x)
}
