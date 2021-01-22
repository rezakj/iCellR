#' iCellR KNN Network
#'
#' This function takes an object of class iCellR and and runs kNet for dimensionality reduction.
#' @param x An object of class iCellR.
#' @param dist.method the distance measure to be used to compute the dissimilarity matrix. This must be one of: "euclidean", "maximum", "mandatattan", "canberra", "binary", "minkowski" or "NULL". By default, distance="euclidean". If the distance is "NULL", the dissimilarity matrix (diss) should be given by the user. If distance is not "NULL", the dissimilarity matrix should be "NULL".
#' @param zoom Adjusting zoom the higher the number the less sensitivity, default = 400.
#' @param data.type Choose between "tsne", "pca", "umap", default = "pca".
#' @param dims PCA dimentions to be use for clustering, default = 1:20.
#' @param joint Run in Combined or joint fashion as in CCCA and CPCA, default = FALSE.
#' @param col.by If return.graph is TRUE the choose the cluster colors.  Choose between "clusters", "conditions".
#' @param return.graph return igraph object, default = FALSE.
#' @param my.seed seed number, default = 1.
#' @param dim.redux Choose between "tsne", "pca", "umap" to unpack the nodes, default = "umap".
#' @param do.redux Perform dim reudx for unpaking the nodes, default = TRUE.
#' @param run.iclust Perform clustering as well (nor recomanded), default = FALSE.
#' @param add.3d Add 3D KNetL as well, default = FALSE.
#' @param layout.2d Choose your 2D layout, default = "layout_nicely".
#' @param layout.3d Choose your 3D layout, default = "layout_with_fr".
#' @return An object of class iCellR.
#' @import progress
#' @importFrom igraph graph.data.frame cluster_louvain modularity membership layout_nicely layout_with_fr V
#' @export
run.knetl <- function (x = NULL,
                    dist.method = "euclidean",
                    zoom = 300,
                    data.type = "pca",
                    dims = 1:20,
                    joint = FALSE,
                    col.by = "clusters",
                    my.seed = 1,
                    layout.2d = "layout_nicely",
                    layout.3d = "layout_with_fr",
                    add.3d = FALSE,
                    dim.redux = "umap",
                    do.redux = TRUE,
                    run.iclust = FALSE,
                    return.graph = FALSE) {
  #
  start_time1 <- Sys.time()
  #  cluster
  if(data.type == "pca") {
    DATA = (t(x@pca.data))[dims, ]
    message(paste("Getting PCA data"))
  }
  if(data.type == "umap") {
    DATA <- t(x@umap.data[1:2])
    message(paste("Getting UMAP data"))
  }
  if(data.type == "tsne") {
    DATA <- t(x@tsne.data)
    message(paste("Getting tSNE data"))
  }
  #####
  message(paste("   Calculating", dist.method,"distance ..."))
  My.distances = as.matrix(dist(t(DATA),dist.method))
  #####
  #####
  k = zoom
  ncells=dim(DATA)[2]
  cell.num = k
  #####
  ### time
  pb <- progress_bar$new(total = ncells,
                         format = "[:bar] :current/:total (:percent) :elapsedfull eta: :eta",
                         clear = FALSE, width= 60)
  if (joint == FALSE) {
    message(paste("    Finding",cell.num, "neighboring cells per cell ..."))
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
  }
###################
  if (joint == TRUE) {
    message(paste("    Finding",cell.num, "jointly neighboring cells per cell ..."))
    MyDF <- as.data.frame(row.names(My.distances))
    colnames(MyDF) <- "IDs"
    KNN1 = lapply(1:ncells, function(findKNN){
      ###### time
      pb$tick()
      #### end time
      ######### loop
      for(i in conditions){
        ha <- paste("^",i,"_",sep="")
        #      CellOrd <- colnames(My.distances)[(GETord(My.distances[,findKNN]))]
        CellOrd <- colnames(My.distances)[(order(My.distances[,findKNN]))]
        CellsId <- subset(CellOrd, grepl(ha, CellOrd))[1:cell.num]
        MyOrd <- as.numeric(rownames(subset(MyDF,MyDF$IDs %in% CellsId)))
        MyDist <- My.distances[MyOrd]
        MyRoot <- rep(findKNN,cell.num)
        data <- cbind(MyRoot,MyOrd,MyDist)
        colnames(data)<- c("from","to","weight")
        MYData <- as.data.frame(data)
        #          CellsId <- grep(ha,CellOrd,value=T, invert=F)[1:myNN]
        NameCol=paste("MySet",i,sep="_")
        eval(call("<-", as.name(NameCol), MYData))
      }
      filenames <- ls(pattern="MySet_")
       TheData <- do.call('rbind', mget(filenames))
       rownames(TheData) <- NULL
       TheData
    })
  }
##############
    data <- do.call("rbind", KNN1)
  #####
  data <- do.call("rbind", KNN1)
    message(paste("    Generating graph from root to neighboring cells ..."))
  g <- graph.data.frame(data, directed=FALSE)
#########
  if (return.graph == FALSE){
    message("Generating 2D Layouts ...")
    set.seed(my.seed)
    data2 <- get(layout.2d)(g)
    colnames(data2) <- c("V1","V2")
    row.names(data2) <- colnames(DATA)
#    data2 <- as.data.frame(data2)
    data2 <- as.data.frame(scale(data2))
    if (add.3d == TRUE) {
      message("Generating 3D Layouts ...")
      data3 <- get(layout.3d)(g, dim =3)
      colnames(data3) <- c("V1","V2","V3")
      data3 <- as.data.frame(scale(data3))
      row.names(data3) <- colnames(DATA)
    }
    ########
    y <- x
    y@pca.data <- data2
  #######
    if (do.redux == TRUE) {
      if (dim.redux == "tsne") {
        message("Running Dimensionality Reduction (tSNE) ...")
        y <- run.pc.tsne(y, dims = 1:2, add.3d = FALSE)
        data2 <- y@tsne.data
      }
      if (dim.redux == "umap") {
        message("Running Dimensionality Reduction (UMAP) ...")
        y <- run.umap(y, dims = 1:2)
        data2 <- y@umap.data
      }
      if (dim.redux == "pca") {
        message("Running Dimensionality Reduction (PCA) ...")
        y@main.data <- as.data.frame(t(data2))
        y <- run.pca(y,scale.data = FALSE)
        data2 <- y@pca.data
      }
    }
    ##################################
    #  png('network_1.png',width = 30, height = 30, units = 'in', res = 300)
    #######
    end_time1 <- Sys.time()
    Time = difftime(end_time1,start_time1,units = "mins")
    Time = round(as.numeric(Time),digits = 2)
    message(paste("Total time",Time,"mins"))
    message(paste("All done!"))
    attributes(x)$knetl.data <- data2
    if (add.3d == TRUE) {
      attributes(x)$knetl.data.3d <- data3
    }
    if (run.iclust == TRUE) {
      message("Running clustering ...")
      y <- iclust(y, sensitivity = k, dims = 1:2)
      x@best.clust <- y@best.clust
    }
    return(x)
  }
  if (return.graph == TRUE){
    ### get clusters
    if (col.by == "clusters") {
      n = max(x@best.clust$clusters)# +1
      getCols <- x@best.clust$clusters
    }
    # get conditions
    if (col.by == "conditions") {
      ha <- colnames(My.distances)
      ha <- data.frame(do.call('rbind', strsplit(as.character(ha),'_',fixed=TRUE)))[1]
      da <- (as.character(as.matrix(ha)))
      getCols <- as.numeric(factor(da))
      conditions <- unique(as.character(as.matrix(ha)))
      n = length(conditions)
    }
    ######## get colors like ggplot
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    ### get colors
    colbar = gg_color_hue(n)
    ########
    V(g)$color <- colbar[getCols]
    end_time1 <- Sys.time()
    Time = difftime(end_time1,start_time1,units = "mins")
    Time = round(as.numeric(Time),digits = 2)
    message(paste("Total time",Time,"mins"))
    message(paste("All done!"))
    return(g)
  }
}
