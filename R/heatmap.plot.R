#' Create heatmaps for genes in clusters or conditions.
#'
#' This function takes an object of class iCellR and genes and provides a heatmap.
#' @param x A data frame containing gene counts for cells.
#' @param gene A set of gene names to be heatmapped.
#' @param cluster.by Choose from "clusters" or "conditions", defult = "clusters".
#' @param cluster.rows If set to FALSE the genes would not be clustered, defult = TRUE.
#' @param scale Choose from "row" or "column", defult = "row".
#' @param heat.colors Colors for heatmap, defult = c("blue" ,"white", "red").
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' MyGenes <- c("SOD1","CD7")
#' MyGenes <- top.markers(marker.genes, topde = 10, min.base.mean = 0.8)
#' heatmap.plot (my.obj, gene = MyGenes, cluster.by = "clusters", cluster.rows = T)
#' }
#' @import pheatmap
#' @export
heatmap.plot <- function (x = NULL,
                          gene = "NULL",
                          cluster.by = "clusters",
                          cluster.rows = F,
                          heat.colors = c("blue","white", "red")) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ## get main data
  DATAmain <- x@main.data
  AllGenes = row.names(DATAmain)
  absent = which((gene %in% AllGenes) == F)
  absentgenes = gene[absent]
  if(length(absentgenes) != 0)
  {
    absentgenes = paste(absentgenes, collapse=",")
    ToPrint <- paste(absentgenes, "not available in your data.
To see the gene names issue this command: row.names(YOURobject@main.data)", sep=" ")
    stop(print(ToPrint))
  }
  ##### get cluster data
      DATA <- x@best.clust
### get main data
  sub.data <- subset(DATAmain,rownames(DATAmain) %in% gene)
  data.t <- t(sub.data)
  data.expr <- as.data.frame(data.t)
## clusters
  if (cluster.by == "clusters") {
  clusters = DATA
# merge
  mrgd <- merge(clusters, data.expr, by="row.names")
  row.names(mrgd) <- mrgd$Row.names
  mrgd <- mrgd[,-1]
  mrgd <- (mrgd[order(as.numeric(as.character(mrgd$clusters)), decreasing = F),])
  SideCol <- mrgd[1]
  SideCol$clusters <- sub("^", "cl.",SideCol$clusters)
  data <- mrgd[,-1]
  data <- as.matrix(data)
  data <- log2(data + 1)
#  data <- t(data)
#  data <- t(scale(data, center = scale.center, scale = TRUE))
  }
  # conditions
  if (cluster.by == "conditions") {
    SideCol <- data.frame(do.call('rbind', strsplit(as.character(rownames(data.expr)),'_',fixed=TRUE)))[1]
    conditions <- as.character(as.matrix(SideCol))
    MyRows = rownames(data.expr)
    conditions <- as.data.frame(cbind(MyRows,conditions = conditions))
    row.names(conditions) <- conditions$MyRows
    conditions <- conditions[2]
    # merge
    mrgd <- merge(conditions, data.expr, by="row.names")
    row.names(mrgd) <- mrgd$Row.names
    mrgd <- mrgd[,-1]
    mrgd <- (mrgd[order(mrgd$conditions, decreasing = F),])
    SideCol <- mrgd[1]
    data <- mrgd[,-1]
    data <- data.matrix(data)
    data <- log2(data + 1)
  }
# Heat map
  mycol <- colorRampPalette(heat.colors)(n = 100)
# fix order
  ClustersOrder <- c(1:length(gene))
  ClustersOrder <- log2(ClustersOrder +1)
#  ClustersOrder <- as.numeric(paste(ClustersOrder,".05", sep=""))
  DATA <- as.data.frame(cbind(gene,ClustersOrder))
  rownames(DATA) <- DATA$gene
  DATA <- (DATA)[2]
  data <- as.data.frame(t(data))
  mrgd <- merge(DATA,data, by="row.names")
  row.names(mrgd) <- mrgd$Row.names
  mrgd <- mrgd[,-1]
#  mrgd <- as.data.frame(as.matrix(mrgd))
  mrgd <- mrgd[order(mrgd$ClustersOrder, decreasing = F),]
  mrgd <- mrgd[,-1]
  data <- t(mrgd)
###########
############ ggplot 2
  # fix order
#  data <- as.data.frame(data)
#  Merged <- merge(SideCol, data, by="row.names")
#  Merged <- Merged[order(nchar(Merged$clusters),Merged$clusters, decreasing = F),]
#  row.names(Merged) <- Merged$Row.names
#  Merged <- Merged[,-1]
#  SideCol <- Merged[1]
#  Merged <- Merged[,-1]
#  data <- as.matrix(Merged)
####### heatmap.2(data, col = mycol,
#  key=TRUE, symm=F,symkey=F,
#  trace="none",density.info="none",
#  Rowv=F, Colv=F,scale="row",
#  cexRow=0.6,cexCol=0.8,margins=c(11,7))
  ###
###
#  row.names(Merged) <- Merged$Row.names
#  Merged <- Merged[,-1]
#  clust.data <- Merged$clusters
#  Merged <- Merged[,-1]
#  Merged <- as.data.frame(t(Merged))
#  Merged <- cbind(gene = rownames(Merged), Merged)
#  Merged <- melt(Merged, id.vars = "gene")
#
#  ggplot(Merged , aes(variable, gene)) +
#    geom_tile(aes(fill = value), color = "white") +
#    scale_fill_gradient(low = "white", high = "steelblue") +
#    ylab("genes ") +
#    xlab("cells") +
#    theme(legend.title = element_text(size = 10),
#          legend.text = element_text(size = 12),
#          plot.title = element_text(size=16),
#          axis.title=element_text(size=14,face="bold"),
#          axis.text.x = element_text(angle = 90, hjust = 1)) +
#    labs(fill = "Expression")
#
# return
  return(pheatmap(t(data),
                  col = mycol,
                  show_colnames = F,
                  cluster_rows = cluster.rows,
                  cluster_cols = F,
                  annotation_col = SideCol,
                  scale = "column"))
}
