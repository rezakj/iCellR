#' iCellR KNN Network
#'
#' This function takes an object of class iCellR and and runs kNet for dimensionality reduction.
#' @param x An object of class iCellR.
#' @param dist.method the distance measure to be used to compute the dissimilarity matrix. This must be one of: "euclidean", "maximum", "mandatattan", "canberra", "binary", "minkowski" or "NULL". By default, distance="euclidean". If the distance is "NULL", the dissimilarity matrix (diss) should be given by the user. If distance is not "NULL", the dissimilarity matrix should be "NULL".
#' @param k KNN the higher the number the less sensitivity, default = 5.
#' @param data.type Choose between "tsne", "pca", "umap", default = "pca". We highly recommend PCA.
#' @param dims PCA dimentions to be use for clustering, default = 1:20.
#' @param my.seed seed number, default = 1.
#' @param conds.to.plot Choose the conditions you want to see in the plot, default = NULL (all conditions).
#' @param my.layout Choose a layout, default = "layout_with_fr".
#' @param node.size Size of the nodes, , default = 10.
#' @param cluster.membership Calculate memberships based on distance.
#' @param abstract Draw all the cells or clusters, , default = TRUE.
#' @param node.colors Color of the nodes, default = random colors.
#' @param edge.color Solor of the edges, default = "gray".
#' @param interactive If set to TRUE an interactive HTML file will be created, default = TRUE.
#' @param out.name If "interactive" is set to TRUE, the out put name for HTML, default = "Abstract.KNetL".
#' @return A plot.
#' @import progress
#' @importFrom igraph graph.data.frame E as.undirected cluster_fast_greedy get.edgelist simplify cluster_louvain modularity membership layout_nicely layout_with_fr V
#' @export
pseudotime.knetl <- function (x = NULL,
                       dist.method = "euclidean",
                       k = 5,
                       abstract = TRUE,
                       data.type = "pca",
                       dims = 1:20,
                       conds.to.plot = NULL,
                       my.layout = "layout_with_fr",
                       node.size = 10,
                       cluster.membership = FALSE,
                       interactive = TRUE,
                       node.colors = NULL,
                       edge.color = "gray",
                       out.name = "Pseudotime.Abstract.KNetL",
                       my.seed = 1
                       ) {
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
  MYtitle = "Pseudotime Abstract KNetL (PAK)"
  MyClusters <- x@best.clust
  if (!is.null(conds.to.plot)) {
    Conditions <- data.frame(do.call('rbind', strsplit(as.character(colnames(DATA)),'_',fixed=TRUE)))[1]
    Conditions <- as.character(as.matrix(Conditions))
    Conditions <- as.data.frame(cbind(colnames(DATA),Conditions))
    colnames(Conditions) <- c("row","conds")
    Conditions <- subset(Conditions, Conditions$conds %in% conds.to.plot)
    Conditions <- as.character(Conditions$row)
    DATA <- as.data.frame(t(DATA))
    DATA <- subset(DATA, rownames(DATA) %in% Conditions)
    DATA <- t(DATA)
    MYtitle = paste(MYtitle," (",conds.to.plot,")", sep="")
    MyClusters <- subset(MyClusters, rownames(MyClusters) %in% Conditions)
  }
  message(paste("   Calculating", dist.method,"distance ..."))
  My.distances = as.matrix(dist(t(DATA),dist.method))
  #####
  #####
  ncells=dim(DATA)[2]
  cell.num = k
  #####  ### time
  pb <- progress_bar$new(total = ncells,
                         format = "[:bar] :current/:total (:percent) :elapsedfull eta: :eta",
                         clear = FALSE, width= 60)
    message(paste("    Finding",cell.num, "neighboring cells per cell ..."))
    KNN1 = lapply(1:ncells, function(findKNN){
      pb$tick()
      #  order(My.distances[,findKNN])[2:cell.num]})
      #  MyOrd <- GETord(My.distances[,findKNN])[2:cell.num]
      MyOrd <- order(My.distances[,findKNN])[2:cell.num]
      MyOrd.IDs <- row.names(My.distances)[MyOrd]
      MyOrd.IDs.clust <- as.numeric(as.matrix(subset(MyClusters, row.names(MyClusters) %in% MyOrd.IDs)))
      ####
      MyDist <- My.distances[MyOrd]
      ####
      MyRoot <- rep(findKNN,cell.num-1)
      MYRoot.ID <- colnames(My.distances)[findKNN]
      MYRoot.IDs <- rep(MYRoot.ID,length(MyRoot))
      MYRoot.clust <- as.numeric(as.matrix(subset(MyClusters, row.names(MyClusters) %in% MYRoot.ID)))
      MYRoot.clusts <- rep(MYRoot.clust,length(MyRoot))
      ####
      data <- cbind(MyRoot,MyOrd,MyDist,MYRoot.IDs,MyOrd.IDs,MYRoot.clusts,MyOrd.IDs.clust)
      colnames(data)<- c("from","to","weight","from.id","to.id","from.clust","to.clust")
      data <- as.data.frame(data)
  })
  ########
    if (abstract == TRUE) {
      message("    Generating abstract graph ...")
      data <- do.call("rbind", KNN1)
      df <- as.data.frame(table(data[6:7]))
      colnames(df) <- c("from","to","weight")
      #    df <- subset(df, df$weight != 0)
      df <- subset(df, df$weight > 5)
      gg <- graph.data.frame(df, directed=FALSE)
      message("   Simplifying graph ...")
      gg <- simplify(gg)
    }
    #######
    if (abstract == FALSE) {
      message("    Generating graph ...")
      data <- do.call("rbind", KNN1)
      df <- as.data.frame(data[1:3])
      colnames(df) <- c("from","to","weight")
      gg <- graph.data.frame(df, directed=FALSE)
    }
    ### make colors
####################
######
   if (is.null(node.colors)) {
     gg_color_hue <- function(n) {
       hues = seq(15, 375, length = n + 1)
       hcl(h = hues, l = 65, c = 100)[1:n]
     }
     n = max(MyClusters$clusters)# +1
     colbar = gg_color_hue(n)
   } else {
     colbar = node.colors
   }
    if (abstract == FALSE) {
      colbar <- colbar[MyClusters$clusters]
    }
###########
#    V(gg)$color <- colbar
    # barplot(c(1:11), col=colbar)
    ####
    MYedgeWith <- log2(as.numeric(E(gg)$weight) + 1)
    # edge.betweenness(gg)
#    ceb <- cluster_edge_betweenness(gg)
#    ceb <- cluster_label_prop(gg)
 #   plot(ceb,gg,layout=layout_with_fr)
## https://briatte.github.io/ggnet/
   message("Drawing layout and generating plot ...")
   if (interactive == TRUE) {
     ####### plotly
     graph = gg
     set.seed(my.seed)
     L <- get(my.layout)(graph)
     vs <- V(graph)
     es <- as.data.frame(get.edgelist(graph))
     Ne <- length(es[1]$V1)
     Xn <- L[,1]
     Yn <- L[,2]
     #### colors
     #   MyCols <- cbind(names(vs),c(1:length(names(vs))), as.character(colbar))
     #   MyCols <- MyCols[order(as.numeric(MyCols[,1]), decreasing = F),]
     #   MyCols <- MyCols[,3]
     #   barplot(c(1:11), col=colbar)
     #   barplot(c(1:11), col=MyCols)
     #####
     edge_shapes <- list()
     for(i in 1:Ne) {
       v0 <- es[i,]$V1
       v1 <- es[i,]$V2
       edge_shape = list(
         type = "line",
         line = list(color = rep(edge.color),
                     width = MYedgeWith[i]),
         x0 = Xn[match(v0,names(vs))],
         y0 = Yn[match(v0,names(vs))],
         x1 = Xn[match(v1,names(vs))],
         y1 = Yn[match(v1,names(vs))],
         opacity = 0.5,
         layer="below"
       )
       edge_shapes[[i]] <- edge_shape}
     #####3
     if (abstract == FALSE) {
       network <- plot_ly(type = "scatter",
                          x = Xn,
                          y = Yn,
                          mode = "markers+text",
                          opacity = 1,
                          marker = list(size = (node.size/2), color = colbar),
                          )
       network <- layout(network,shapes = edge_shapes,
                         xaxis = list(title = "", showgrid = FALSE,
                                      showticklabels = FALSE, zeroline = FALSE),
                         yaxis = list(title = "", showgrid = FALSE,
                                      showticklabels = FALSE, zeroline = FALSE)) %>%
         layout(plot_bgcolor = "white") %>%
         layout(paper_bgcolor = "white") %>%
         layout(title = MYtitle)
     }
     if (abstract == TRUE) {
       network <- plot_ly(type = "scatter",
                          x = Xn,
                          y = Yn,
                          mode = "markers+text",
                          opacity = 1,
                          marker = list(size = (node.size * 2), color = colbar),
                          text = names(vs),)
       #######
       network <- layout(network,shapes = edge_shapes,
                         xaxis = list(title = "", showgrid = FALSE,
                                      showticklabels = FALSE, zeroline = FALSE),
                         yaxis = list(title = "", showgrid = FALSE,
                                      showticklabels = FALSE, zeroline = FALSE)) %>%
         layout(plot_bgcolor = "white") %>%
         layout(paper_bgcolor = "white") %>%
         layout(title = MYtitle)
     }
     ######
     OUT.PUT <- paste(out.name, ".html", sep="")
     htmlwidgets::saveWidget(ggplotly(network), OUT.PUT)
   } else {
   ###############
   rm(.Random.seed, envir=.GlobalEnv)
   set.seed(my.seed)
   if (cluster.membership == TRUE) {
     ceb <- cluster_fast_greedy(as.undirected(gg))
     return(
       plot(ceb,gg, layout=get(my.layout),
            vertex.size=node.size,
            vertex.color=colbar,
            edge.width=MYedgeWith,
            main = MYtitle,
            edge.color=edge.color))
   }
#######
   if (cluster.membership == FALSE) {
     if (abstract == FALSE) {
       return(
       plot(gg, layout=layout_nicely,
            vertex.size=(node.size/2),
            vertex.label=NA,
            margin=0,
            vertex.color=colbar))
     }
     return(
       plot(gg, layout=get(my.layout),
            vertex.size=node.size,
            vertex.color=colbar,
            edge.width=MYedgeWith,
            main = MYtitle,
            edge.color=edge.color))
    }
  }
}

