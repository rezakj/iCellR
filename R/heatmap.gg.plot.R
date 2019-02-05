#' Create heatmaps for genes in clusters or conditions.
#'
#' This function takes an object of class iCellR and genes and provides a heatmap.
#' @param x A data frame containing gene counts for cells.
#' @param gene A set of gene names to be heatmapped.
#' @param cluster.by Choose from "clusters" or "conditions", defult = "clusters".
#' @param heat.colors Colors for heatmap, defult = c("blue" ,"white", "red").
#' @param interactive If TRUE an html intractive file will be made, defult = TRUE.
#' @param out.name Output name for html file if interactive = TRUE defult = "plot".
#' @param no.key If you want a color legend key defult = FALSE.
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' heatmap.gg.plot(my.obj, gene = MyGenes, interactive = T, out.name = "plot", cluster.by = "clusters")
#' }
#' @import pheatmap
#' @importFrom reshape melt
#' @export
heatmap.gg.plot <- function (x = NULL,
                          gene = "NULL",
                          cluster.by = "clusters",
                          min.scale = -2.5,
                          max.scale = 2.5,
                          interactive = T,
                          cex.col = 10,
                          cex.row = 10,
                          no.key = FALSE,
                          out.name = "plot",
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
  MYord <- cbind(Row = rownames(DATA), DATA)
  ## order by cluster
  if (cluster.by == "clusters") {
#  MYord <- (MYord[order(MYord$Row, decreasing = F),])
  MYord <- (MYord[order(MYord$clusters, decreasing = F),])
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
  if (cluster.by == "conditions") {
    MYord <- (MYord[order(MYord$Row, decreasing = F),])
    z = as.data.frame(MYord$Row)
  }
  # order
  MYord <- row.names(MYord)
  ### get main data
  sub.data <- DATAmain[gene, colnames(DATAmain),drop = FALSE]
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
  ############ ggplot 2
  # fix order
  library(reshape)
  data <- melt(data)
#  head(data)
  names(x = data)[names(x = data) == "X1"] <- "cell"
  names(x = data)[names(x = data) == "X2"] <- "gene"
  names(x = data)[names(x = data) == "value"] <- "expression"
  head(data)
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
  data <- cbind(Myord = row.names(data), data)
  mrgd <- merge(data, clusters, by="cell")
  mrgd <- (mrgd[order(as.numeric(as.character(mrgd$Myord)), decreasing = T),])
  data <- mrgd
  data$gene <- factor(data$gene, levels = rev(gene))
  data$cell <- factor(data$cell, levels = MYord)
  ### plot
  heatmap <- ggplot(data, aes(x = cell, y = gene, fill = expression, text=clusters)) + geom_tile() +
    scale_fill_gradient2(low = col.low, mid = col.mid, high = col.high, name = "Expression",
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

  heatmap <- heatmap + facet_wrap( ~ clusters, nrow = 1,scales = "free_x") +
    theme(strip.background = element_rect(fill= NA))

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
  if (no.key == T) {
    heatmap <- heatmap + theme(legend.position = "none")
  }
#
#  htmlwidgets::saveWidget(ggplotly(heatmap), "fix.html")
  # return
  if (interactive == T) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(heatmap), OUT.PUT)
  } else {
    return(heatmap)
  }
}
