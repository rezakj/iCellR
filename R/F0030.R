#' Create heatmaps for genes in clusters or conditions.
#'
#' This function takes an object of class iCellR and genes and provides a heatmap.
#' @param x A data frame containing gene counts for cells.
#' @param gene A set of gene names to be heatmapped.
#' @param cluster.by Choose from "clusters" or "none", default = "clusters".
#' @param conds.to.plot Choose the conditions you want to see in the plot, default = NULL (all conditions).
#' @param heat.colors Colors for heatmap, default = c("blue" ,"white", "red").
#' @param cell.sort If FALSE the cells will not be sorted based on their distance, default = TRUE.
#' @param interactive If TRUE an html interactive file will be made, default = TRUE.
#' @param out.name Output name for html file if interactive = TRUE, default = "plot".
#' @param no.key If you want a color legend key, default = FALSE.
#' @param data.type Choose from "main" and "imputed", default = "main".
#' @param min.scale Set a minimum color scale, default = -2.5.
#' @param max.scale Set a maximum color scale, default = 2.5.
#' @param cex.col Chhose a size, default = 10.
#' @param cex.row Choose a size, default = 10.
#' @return An object of class iCellR
#' @examples
#' marker.genes <- findMarkers(demo.obj,fold.change = 2,padjval = 0.1,uniq = TRUE)
#'
#' MyGenes <- top.markers(marker.genes, topde = 10, min.base.mean = 0.8)
#'
#' heatmap.gg.plot(demo.obj,
#'                gene = MyGenes,
#'                out.name = "plot",
#'                cluster.by = "clusters",
#'                interactive = FALSE)
#' @import pheatmap
#' @importFrom reshape melt
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly ggplotly layout plot_ly
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @importFrom methods new
#' @importFrom stats aggregate as.dendrogram cor cor.test dist hclust p.adjust prcomp quantile sd t.test
#' @importFrom utils capture.output packageVersion read.table write.table
#' @importFrom graphics legend par plot
#' @importFrom ggplot2 ggplot geom_segment geom_violin guide_colorbar guide_legend guides scale_color_discrete scale_colour_gradient scale_fill_gradient2 scale_x_continuous scale_y_continuous scale_y_discrete stat_summary coord_polar element_rect element_text element_blank facet_wrap scale_color_manual geom_hline geom_jitter geom_vline ylab xlab ggtitle theme_bw aes theme geom_bar geom_point geom_boxplot geom_errorbar position_dodge geom_tile geom_density geom_line
#' @export
heatmap.gg.plot <- function (x = NULL,
                          gene = "NULL",
                          cell.sort = FALSE,
                          data.type = "main",
                          cluster.by = "clusters",
                          conds.to.plot = NULL,
                          min.scale = -2.5,
                          max.scale = 2.5,
                          interactive = TRUE,
                          cex.col = 10,
                          cex.row = 10,
                          no.key = FALSE,
                          out.name = "plot",
                          heat.colors = c("blue","white", "red")) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ## get main data
  if (data.type == "main") {
    DATAmain <- x@main.data
  }
  if (data.type == "imputed") {
    DATAmain <- x@imputed.data
  }
  AllGenes = row.names(DATAmain)
  absent = which((gene %in% AllGenes) == FALSE)
  absentgenes = gene[absent]
  #######
  DATAmain <- DATAmain[ rowSums(DATAmain) > 0, ]
  AllGenes2 = row.names(DATAmain)
  Gene0 <- setdiff(AllGenes,AllGenes2)
  Gene0List = which((gene %in% Gene0) == TRUE)
  Gene0List = gene[Gene0List]
  NotInHeatmap <- c(absentgenes,Gene0List)
  gene <- setdiff(gene,NotInHeatmap)
#####
  if(length(absentgenes) != 0)
  {
    absentgenes = paste(absentgenes, collapse=",")
    ToPrint <- paste("WARNING:",absentgenes, "not available in your data", sep=" ")
    message(ToPrint)
    Gene0List = paste(Gene0List, collapse=",")
    ToPrint1 <- paste("WARNING:",Gene0List, "not expressed in your data", sep=" ")
    message(ToPrint1)
    presentgenes = paste(gene, collapse=",")
    ToPrint2 <- paste("WARNING:","plotting",presentgenes, sep=" ")
    message(ToPrint2)
  }
  ##### get cluster data
  DATA <- x@best.clust
  MYord <- cbind(Row = rownames(DATA), DATA)
 ############
  ## order by cluster
  if (cluster.by == "clusters") {
#  MYord <- (MYord[order(MYord$Row, decreasing = F),])
  MYord <- (MYord[order(MYord$clusters, decreasing = FALSE),])
  clustOrd <- unique(MYord$clusters)
  z = as.data.frame(MYord$clusters)
  # make break lines
  for(i in 1:length(clustOrd)) {
    NameCol=paste("cluster",i,sep="_")
    myDATA = as.numeric(tail(row.names(subset(z, MYord$clusters == i)),1))
    eval(call("<-", as.name(NameCol), myDATA))
  }
  MyLines <- as.numeric(mget(ls(pattern="cluster_")))
  MyLines <- MyLines[1:length(MyLines)-1]
  }
  ################
  if (cluster.by == "conditions") {
    clustOrd <-  data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
    colnames(clustOrd)="clusters"
    MYord$clusters <- clustOrd
    cond.data <- MYord
    ha <- cond.data$clusters
    colnames(ha) <- "MyConds"
    cond.data$MyConds <- ha$MyConds
    z <- as.data.frame(MYord$clusters)
    clustOrd = unique(MYord$clusters)
    MyLines <- as.numeric(row.names(clustOrd))
  }
  # order
  MYord <- row.names(MYord)
  ### get main data
  sub.data <- DATAmain[gene, colnames(DATAmain),drop = FALSE]
  if (cell.sort == TRUE) {
    counts.pca <- prcomp(t(scale(t(sub.data))), center = FALSE, scale. = FALSE)
    my.data.my.pca = t(counts.pca$rotation)[1:2, ]
    My.distances = as.data.frame(as.matrix(dist(t(my.data.my.pca))))[1]
    colnames(My.distances) <- "MyDist"
    My.distances$MyIDs <- row.names(My.distances)
    My.distances <- row.names(My.distances[order(My.distances$MyDist, decreasing = T),])
  }
##########
  data.t <- t(sub.data)
  data.expr <- as.data.frame(data.t)
  ## clusters
    clusters = DATA
    cl.cell <- paste(DATA$clusters, row.names(DATA), sep = "_")
    data <- as.matrix(data.expr)
    # data <- scale(log2(data + 1)) # this makes better cols but fails with moothing
     data <- scale(data)
    # fix scale
    FixScale <- function (mydata, min, max){
      Mydat <- mydata
      Mydat[Mydat > max] <- max
      Mydat[Mydat < min] <- min
      return(Mydat)
    }
    #
    data <- FixScale(mydata = data, min = min.scale, max = max.scale)
  ###########
  col.low = heat.colors[1]
  col.mid = heat.colors[2]
  col.high = heat.colors[3]
  ############ get the conditions you need
  if (!is.null(conds.to.plot)) {
    MyCondsIndata <-  data.frame(do.call('rbind', strsplit(as.character(rownames(data)),'_',fixed=TRUE)))[1]
    MyCondsIndata <- cbind(MyCondsIndata,rownames(data))
    colnames(MyCondsIndata) <- c("conds","MYrows")
    MyCondsIndata <- subset(MyCondsIndata, MyCondsIndata$conds %in% conds.to.plot)
    MyCondsIndata <- MyCondsIndata$MYrows
#    dim(data)
    data <- subset(data, row.names(data) %in% MyCondsIndata)
#    dim(data)
  }
  ############ ggplot 2
  data <- melt(data)
#  head(data)
  names(x = data)[names(x = data) == "X1"] <- "cell"
  names(x = data)[names(x = data) == "X2"] <- "gene"
  names(x = data)[names(x = data) == "value"] <- "expression"
  # head(data)
  ## Make color breaks based on expression
  CustomPalette <- function (low = "white", high = "red", mid = NULL, k = 50)
  {
    low <- col2rgb(col = low)/255
    high <- col2rgb(col = high)/255
    if (is.null(x = mid)) {
      r <- seq(from = low[1], to = high[1], len = k)
      g <- seq(from = low[2], to = high[2], len = k)
      b <- seq(from = low[3], to = high[3], len = k)
    }
    else {
      k2 <- round(x = k/2)
      mid <- col2rgb(col = mid)/255
      r <- c(seq(from = low[1], to = mid[1], len = k2), seq(from = mid[1],
                                                            to = high[1], len = k2))
      g <- c(seq(from = low[2], to = mid[2], len = k2), seq(from = mid[2],
                                                            to = high[2], len = k2))
      b <- c(seq(from = low[3], to = mid[3], len = k2), seq(from = mid[3],
                                                            to = high[3], len = k2))
    }
    return(rgb(red = r, green = g, blue = b))
  }
#
  My.50.col <- function (...)
  {
    return(CustomPalette(low = "blue", high = "white", mid = "red",
                         ...))
  }
# Final color breaks
  breaks <- seq(from = min(data$expression),
                to = max(data$expression),length = length(x = My.50.col()) + 1)
  # genes
  data$gene <- with(data = data,
                    expr = factor(x = gene,levels = rev(x = unique(x = data$gene))))
#
  # merge data with cluster
  clusters <- cbind(cell = rownames(clusters),clusters)
  if (cluster.by == "conditions") {
    colnames(cond.data) <- c("cell","clusters","MyConds")
    clusters <- cond.data
  }
  data <- cbind(Myord = row.names(data), data)
  mrgd <- merge(data, clusters, by="cell")
  mrgd <- (mrgd[order(as.numeric(as.character(mrgd$Myord)), decreasing = TRUE),])
  data <- mrgd
  data$gene <- factor(data$gene, levels = rev(gene))
  data$cell <- factor(data$cell, levels = MYord)
  if (cell.sort == TRUE) {
    data$cell <- factor(data$cell, levels = My.distances)
  }
  ### plot
  heatmap <- ggplot(data, aes(x = cell, y = gene, fill = expression, text=clusters)) + geom_tile() +
    scale_fill_gradient2(low = col.low, mid = col.mid, high = col.high, name = "",
    guide = guide_colorbar(direction = "vertical", title.position = "left")) +
    scale_y_discrete(position = "right", labels = rev(gene)) +
    theme(axis.line = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), strip.text.x = element_text(size = 15),
          axis.text.y = element_text(size = cex.row), axis.text.x = element_text(size = cex.col),
          axis.title.x = element_blank())
#
  heatmap <- heatmap + theme(axis.title.x = element_blank(),
                             axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                             axis.line = element_blank(), axis.title.y = element_blank(),
                             axis.ticks.y = element_blank())
  #
#  Clust.Ord <- factor(data$clusters, levels = clustOrd)
#  heatmap <- heatmap + scale_x_discrete(labels = Clust.Ord, position = "top")
#  heatmap <- heatmap + scale_x_discrete(position = "top")
#  heatmap <- heatmap + scale_x_discrete(labels = as.character(sort(as.numeric(Clust.Ord))))
#  heatmap <- heatmap + scale_x_discrete(breaks=as.character(sort(as.numeric(Clust.Ord))),
#                        labels=as.character(sort(as.numeric(Clust.Ord))),position = "top")
#  geom_vline(xintercept = top.rank.line, colour = rank.line.col) +
#  geom_hline(yintercept = SDlimit, colour = disp.line.col) +
#  my.lines<-data.frame(x=c(.5,4.5), y=c(5.5,.5), xend=c(4.5,4.5), yend=c(5.5,5.5))
#  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), size=3, inherit.aes=F)

#########
#  if (clustline == T) {
#    heatmap <- heatmap + geom_vline(xintercept = MyLines,
#                                    alpha = line.transparency,
#                                    size = line.size,
#                                    colour = "black")
#  }

  ###
  if (cluster.by == "clusters") {
  heatmap <- heatmap + facet_wrap( ~ clusters, nrow = 1,scales = "free_x") +
    theme(strip.background = element_rect(fill= NA))
  }
  if (cluster.by == "conditions") {
    heatmap <- heatmap + facet_wrap( ~ MyConds, nrow = 1,scales = "free_x") +
      theme(strip.background = element_rect(fill= NA))
  }
  if (cluster.by == "none") {
    heatmap <- heatmap
  }

  # line per group of genes
#  if (geneline == T) {
#    PerClust <- length(gene) / tail(clustOrd,1)
#    PerClust <- (clustOrd * PerClust)
#    PerClust <- PerClust[1:length(PerClust)-1]
#    PerClust <- PerClust + 0.5
#    heatmap <- heatmap + geom_hline(yintercept = PerClust, colour = "black")
#  }
  #
#  panel.spacing <- unit(x = 0.15, units = "lines")
#  heatmap <- heatmap + theme(strip.background = element_blank(),panel.spacing = panel.spacing)
  # rmove key
  if (no.key == TRUE) {
    heatmap <- heatmap + theme(legend.position = "none")
  }
#
#  htmlwidgets::saveWidget(ggplotly(heatmap), "fix.html")
  # return
  if (interactive == TRUE) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(heatmap), OUT.PUT)
  } else {
    return(heatmap)
  }
}
