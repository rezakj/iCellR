#' Pseudotime Tree
#'
#' This function takes an object of class iCellR and marker genes for clusters and performs pseudotime for differentiation or time course analysis.
#' @param x An object of class iCellR.
#' @param clust.method Choose from "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid", default = "complete".
#' @param dist.method Choose from "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski", default = "euclidean".
#' @param clust.names A list of names for clusters.
#' @param marker.genes A list of marker genes for clusters.
#' @param label.offset Space between names and tree, default = 0.5.
#' @param hang Hang, default = 1.
#' @param type Choose from "classic", "jitter", "unrooted", "fan", "cladogram", "radial", default = "classic".
#' @param cex Text size, default = 1.
#' @return An object of class iCellR.
#' @examples
#' marker.genes <- findMarkers(demo.obj,fold.change = 2,padjval = 0.1,uniq = TRUE)
#'
#' MyGenes <- top.markers(marker.genes, topde = 10, min.base.mean = 0.8)
#'
#' pseudotime.tree(demo.obj,
#'                marker.genes = MyGenes,
#'                type = "unrooted",
#'                clust.method = "complete")
#'
#' @import gridExtra
#' @import ggdendro
#' @import ape
#' @importFrom ggplot2 ggplot theme_classic geom_segment geom_violin guide_colorbar guide_legend guides scale_color_discrete scale_colour_gradient scale_fill_gradient2 scale_x_continuous scale_y_continuous scale_y_discrete stat_summary coord_polar element_rect element_text element_blank facet_wrap scale_color_manual geom_hline geom_jitter geom_vline ylab xlab ggtitle theme_bw aes theme geom_bar geom_point geom_boxplot geom_errorbar position_dodge geom_tile geom_density geom_line
#' @export
pseudotime.tree <- function (x = NULL,
                             marker.genes = "NULL",
                             clust.names = "NULL",
                             dist.method = "euclidean",
                             clust.method = "complete",
                             label.offset = 0.5,
                             type = "classic",
                             hang = 1,
                             cex = 1) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  # geth the genes and scale them based on model
  DATA <- x@clust.avg
  row.names(DATA) <- DATA$gene
  DATA <- DATA[,-1]
  if (clust.names[1] != "NULL") {
    colnames(DATA) <- clust.names
  }
  if (marker.genes[1] == "NULL") {
    stop("provide marker genes for clusters (e.g. top 10 for each cluster)")
  }
  MyGenes <- marker.genes
  topGenes <- as.matrix(subset(DATA, rownames(DATA) %in% MyGenes))
  DATA <- dist(scale(t(DATA)), method = dist.method)
  hc <- hclust(DATA, method = clust.method)
##### gitter plot
#  DATA <- x@pca.data[1:5]
#  DATA <- x@tsne.data.3d
  DATA <- x@tsne.data
#  DATA <- x@diffusion.data
  data <- data.matrix(DATA)
  data <- dist(data, method = dist.method)
#  hcgg <- hclust(data, method = clust.method)
#  dhc <- as.dendrogram(hcgg)
  MyPC <- as.data.frame(as.matrix(data))[1]
  colnames(MyPC) <- "distance"
  MyClust <- x@best.clust
  MyGitterData <- cbind(MyPC, MyClust)
  MyOrd <- hc$order
  MyGitterData$clusters <- factor(MyGitterData$clusters, levels = MyOrd)
  dhc <- as.dendrogram(hc)
  ddata <- dendro_data(dhc)
  dend <- ddata$segments
  ##
  P1 <- ggplot(MyGitterData,aes(y=scale(distance),x=as.factor(clusters),col=clusters)) +
    geom_jitter() +
#    geom_violin(trim=FALSE, alpha = 0.5) +
    theme(panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.key = element_rect(fill = "white")) +
    theme(legend.position = "none") + theme(axis.text.x=element_text(angle=90))
  ##
  P2 <- ggplot(dend) +
  geom_segment(aes(x = x,
                     y = y,
                     xend = xend,
                     yend = yend)) +
    theme(panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.key = element_rect(fill = "white")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab("distance")
#####
  if (type == "classic") {
    return(plot(hc, hang = hang, ylab = "Height", xlab = "Clusters", sub=""))
  }
  if (type == "jitter") {
    return(grid.arrange(P2,P1, nrow = 2, heights=c(1,4)))
  }
  if (type != "jitter" || type != "classic") {
    return(plot(as.phylo(hc),
         type = type,
         cex = cex,
         label.offset = label.offset))
  }
}
