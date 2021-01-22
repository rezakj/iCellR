#' iCellR Clustering
#'
#' This function takes an object of class iCellR and finds optimal number of clusters and clusters the data.
#' @param x An object of class iCellR.
#' @param dist.method the distance measure to be used to compute the dissimilarity matrix. This must be one of: "euclidean", "maximum", "mandatattan", "canberra", "binary", "minkowski" or "NULL". By default, distance="euclidean". If the distance is "NULL", the dissimilarity matrix (diss) should be given by the user. If distance is not "NULL", the dissimilarity matrix should be "NULL".
#' @param sensitivity The higher the number the less sensitivity, default = 100.
#' @param data.type Choose between "tsne", "pca", "umap", default = "pca".
#' @param dims PCA dimentions to be use for clustering, default = 1:10.
#' @param return.graph return igraph object, default = FALSE.
#' @return An object of class iCellR.
#' @import progress
#' @importFrom igraph graph.data.frame cluster_louvain modularity membership
#' @export
iclust <- function (x = NULL,
                      dist.method = "euclidean",
                      sensitivity = 100,
                      data.type = "pca",
                      dims = 1:10,
                      return.graph = FALSE) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  #
  start_time1 <- Sys.time()
  #  cluster
 if(data.type == "pca") {
   DATA = (t(x@pca.data))[dims, ]
   message(paste("Clustering based on PCA data"))
 }
 if(data.type == "umap") {
   DATA <- t(x@umap.data[1:2])
   message(paste("Clustering based on UMAP data"))
 }
 if(data.type == "tsne") {
   DATA <- t(x@tsne.data)
   message(paste("Clustering based on tSNE data"))
 }
  if(data.type == "knetl") {
    DATA <- t(x@knetl.data)
    message(paste("Clustering based on KNetL data"))
  }
#####
message(paste("   Calculating", dist.method,"distance ..."))
My.distances = as.matrix(dist(t(DATA),dist.method))
#####
k = sensitivity
ncells=dim(DATA)[2]
cell.num = k
#####
message(paste("    Finding",cell.num, "neighboring cells per cell ..."))
### time
pb <- progress_bar$new(total = ncells,
                       format = "[:bar] :current/:total (:percent) :elapsedfull eta: :eta",
                       clear = FALSE, width= 60)
###
KNN1 = lapply(1:ncells, function(findKNN){
  pb$tick()
#  order(My.distances[,findKNN])[2:cell.num]})
#  MyOrd <- GETord(My.distances[,findKNN])[2:cell.num]
  MyOrd <- order(My.distances[,findKNN])[2:cell.num]
  MyDist <- My.distances[MyOrd]
  MyRoot <- rep(findKNN,cell.num-1)
  data <- cbind(MyRoot,MyOrd,MyDist)
  colnames(data)<- c("from","to","weight")
  data <- as.data.frame(data)
  })
#####
data <- do.call("rbind", KNN1)
#####
### time
#message(paste("    Finding root cells for neighboring cells ..."))
#pb <- progress_bar$new(total = ncells,
#                       format = "[:bar] :current/:total (:percent) :elapsedfull eta: :eta",
#                       clear = FALSE, width= 60)
#KNN2 = lapply(1:ncells, function(findKNN){
#  pb$tick()
#  rep(findKNN,cell.num-1)})
####
### time
#message(paste("    Finding weight for neighboring cells ..."))
#pb <- progress_bar$new(total = ncells,
#                       format = "[:bar] :current/:total (:percent) :elapsedfull eta: :eta",
#                       clear = FALSE, width= 60)
#KNN3 = lapply(1:ncells, function(findKNN){
#  da=order(My.distances[,findKNN])[2:cell.num]
#  pb$tick()
#  My.distances[da]
#})
######### make data frame for igraph
#data = as.data.frame(cbind(unlist(KNN2),unlist(KNN1),unlist(KNN3)))
#colnames(data)<- c("from","to","weight")
### run igraph
message(paste("    Generating graph from root to neighboring cells ..."))
g <- graph.data.frame(data, directed=FALSE)
message(paste("    Clustering ..."))
community <- cluster_louvain(g)
#comps <- membership(community)
#colbar <- rainbow(max(comps)+1)
#V(g)$color <- colbar[comps+1]
clusters = as.numeric(membership(community))
############ export
clustInfo <- as.data.frame(cbind(colnames(DATA),clusters))
row.names(clustInfo) <- clustInfo$V1
clustInfo <- clustInfo[-1]
clustInfo$clusters <- as.numeric(clustInfo$clusters)
message(paste("Found",max(unique(clustInfo$clusters)),"clusters"))
end_time1 <- Sys.time()
Time = difftime(end_time1,start_time1,units = "mins")
Time = round(as.numeric(Time),digits = 2)
message(paste("Total time",Time,"mins"))
message(paste("All done!"))
if (return.graph == FALSE){
attributes(x)$best.clust <- clustInfo
return(x)
}
if (return.graph == TRUE){
  return(g)
 }
}



