#' Make scatter, box and bar plots for genes
#'
#' This function takes an object of class iCellR and provides plots for genes.
#' @param x An object of class iCellR.
#' @param gene A gene name to be plotted.
#' @param box.to.test A cluster number so that all the boxes in the box plot would be compared to. If set to "0" the cluster with the highest avrage would be choosen, default = 0.
#' @param box.pval Choose from "sig.values" and "sig.signs". If set to "sig.signs" p values would be replaced with signs ("na", "*", "**", "***"), default = "sig.signs".
#' @param plot.data.type Choose from "tsne" and "pca", default = "tsne".
#' @param clust.dim 2 for 2D plots and 3 for 3D plots, default = 2.
#' @param col.by Choose from "clusters" and "conditions", default = "clusters".
#' @param cond.shape If TRUE the conditions will be shown in shapes.
#' @param plot.type Choose from "scatterplot", "boxplot" and "barplot", default = "scatterplot".
#' @param cell.size A number for the size of the points in the plot, default = 1.
#' @param cell.colors Colors for heat mapping the points in "scatterplot", default = c("gray","red").
#' @param box.cell.col A color for the points in the box plot, default = "black".
#' @param box.color A color for the boxes in the "boxplot", default = "red".
#' @param box.line.col A color for the lines around the "boxplot", default = "green".
#' @param back.col A color for the plot background, default = "black".
#' @param cell.transparency Color transparency for points in "scatterplot" and "boxplot", default = 0.5.
#' @param interactive If set to TRUE an interactive HTML file will be created, default = TRUE.
#' @param out.name If "interactive" is set to TRUE, the out put name for HTML, default = "plot".
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' cluster.plot(my.obj,
#'             cell.size = 1,
#'             plot.type = "tsne",
#'             cell.color = "black",
#'             back.col = "white",
#'             col.by = "clusters",
#'             cell.transparency = 0.5,
#'             clust.dim = 3,
#'             interactive = T,
#'             density = F,
#'             out.name = "tSNE_3D_clusters")
#'
#' cluster.plot(my.obj, cell.size = 1,
#'             plot.type = "tsne",
#'             col.by = "clusters",
#'             cell.transparency = 0.5,
#'             clust.dim = 2,
#'             interactive = T,
#'             density = F,
#'             out.name = "tSNE_2D_clusters")
#'
#'cluster.plot(my.obj,
#'            cell.size = 2,
#'            plot.type = "tsne",
#'            clust.dim = 2,
#'            interactive = F)
#'
#'cluster.plot(my.obj,
#'            cell.size = 1,
#'            plot.type = "tsne",
#'            col.by = "clusters",
#'            clust.dim = 3,
#'            interactive = F,
#'            angle = 45)
#'
#'
#'cluster.plot(my.obj,
#'           cell.size = 1,
#'           plot.type = "pca",
#'           cell.color = "black",
#'           back.col = "white",
#'           col.by = "conditions",
#'           cell.transparency = 0.5,
#'           clust.dim = 3,
#'           interactive = T,
#'           density = F,
#'           out.name = "PCA_3D_conditions")
#' }
#' @import ggpubr
#' @export
gene.plot <- function (x = NULL,
                       gene = "NULL",
                       cond.shape = F,
                       data.type = "main",
                       box.to.test = 0,
                       box.pval = "sig.signs",
                       plot.data.type = "tsne",
                       min.scale = -1,
                       max.scale = 1,
                       clust.dim = 2,
                       col.by = "clusters",
                       plot.type = "scatterplot",
                       cell.size = 1,
                       cell.colors = c("gray","red"),
                       box.cell.col = "black",
                       box.color = "red",
                       box.line.col = "green",
                       back.col = "white",
                       cell.transparency = 0.5,
                       interactive = TRUE,
                       out.name = "plot") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  if (gene == "NULL") {
    stop("There is no gene name provided. Please provide a gene name")
  }
  if (length(gene) != 1) {
    stop("currently you can only plot one gene at a time")
  }
  ## get main data
  if (data.type == "main") {
  DATAmain <- x@main.data
  }
  if (data.type == "imputed") {
    DATAmain <- x@imputed.data
  }
  AllGenes = row.names(DATAmain)
  gene.availability = gene %in% AllGenes
  if(gene.availability != TRUE)
  {
    stop("Your gene name is not in the main data. To see the gene names issue this command:
         row.names(YOURobject@main.data)")
  }
  ##### get cluster data
  # 2 dimentions
  if (clust.dim == 2) {
    if (plot.data.type == "tsne") {
      MyTitle = "tSNE Plot"
      DATA <- x@tsne.data
    }
    if (plot.data.type == "pca") {
      MyTitle = "PCA Plot"
      DATA <- x@pca.data
    }
    if (plot.data.type == "umap") {
      MyTitle = "UMAP Plot"
      DATA <- x@umap.data
    }
    if (plot.data.type == "diffusion") {
      MyTitle = "Diffusion Map Plot"
      DATA <- x@diffusion.data
    }
  }
  # 3 dimentions
  if (clust.dim == 3) {
    if (plot.data.type == "tsne") {
      MyTitle = "3D tSNE Plot"
      DATA <- x@tsne.data.3d
    }
    if (plot.data.type == "pca") {
      MyTitle = "3D PCA Plot"
      DATA <- x@pca.data
    }
    if (plot.data.type == "dst") {
      MyTitle = "3D DST Plot"
      DATA <- x@diff.st.data
    }
    if (plot.data.type == "diffusion") {
      MyTitle = "Diffusion Map Plot"
      DATA <- x@diffusion.data
    }
  }
  # conditions
  if (col.by == "conditions") {
    col.legend <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
    col.legend.box <- factor(as.matrix(col.legend))
  }
  # clusters
  if (col.by == "clusters") {
    if (is.null(x@best.clust)) {
      stop("Clusters are not assigend yet")
    } else {
      col.legend.box <- x@best.clust
      col.legend.box$clusters <- sub("^", "cl.",col.legend.box$clusters)
      col.legend.box <- factor(col.legend.box$clusters)
    }
  }
###### make binary
  # get the gene from the main data
  sub.data <- subset(DATAmain,rownames(DATAmain) == gene)
  data.t <- t(sub.data)
  data.expr <- as.data.frame(data.t)
  data.binary <- data.t > 0
  data.binary <- as.data.frame(data.binary)
  col.legend.bin = data.binary[[gene]]
#### make heamap
  col.legend = log2(data.expr + 1)
  col.legend <- as.numeric(as.matrix(col.legend))
# fix scale
  FixScale <- function (mydata, min, max){
    Mydat <- mydata
    Mydat[Mydat > max] <- max
    Mydat[Mydat < min] <- min
    return(Mydat)
  }
  col.legend.scaled <- FixScale(mydata = col.legend, min = min.scale, max = max.scale)
# fix color for ADTs
#      if ( length(grep("^ADT_", gene, value = T)) == 1) {
#        Lo3=(quantile(col.legend,0.05)) # 1
#        Lo2=(quantile(col.legend,0.10)) # 2
#        Lo1=(quantile(col.legend,0.25)) # 3
#        MID=(quantile(col.legend,0.50)) # 4
#        Up1=(quantile(col.legend,0.75)) # 5
#        Up2=(quantile(col.legend,0.90)) # 6
#        Up3=(quantile(col.legend,0.95)) # 7
#        col.legend <- replace(col.legend, col.legend > Up3, Up3)
#        col.legend <- replace(col.legend, col.legend > Up2 & col.legend < Up3 ,Up2)
#        col.legend <- replace(col.legend, col.legend > Up1 & col.legend < Up2 ,Up1)
#        col.legend <- replace(col.legend, col.legend > MID & col.legend < Up1 ,MID)
#        col.legend <- replace(col.legend, col.legend > Lo1 & col.legend < MID ,Lo1)
#        col.legend <- replace(col.legend, col.legend < Lo2 ,Lo1)
#      }
#      if ( length(grep("^ADT_", gene, value = T)) == 0) {
#        col.legend = col.legend
#      }
###
  if (plot.type == "scatterplot") {
  # plot 2d
    Conditions = col.legend.box
  if (clust.dim == 2) {
    if (interactive == F) {
      if (cond.shape == F) {
        myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                   text = row.names(DATA), color = col.legend.scaled)) +
          geom_point(size = cell.size, alpha = cell.transparency) +
          scale_colour_gradient(low = cell.colors[1], high = cell.colors[2], name="") +
          xlab("Dim1") +
          ylab("Dim2") +
          ggtitle(paste(MyTitle,"for (",gene,")")) +
          theme(panel.background = element_rect(fill = back.col, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.key = element_rect(fill = back.col))
      }
        if (cond.shape == T) {
          myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                     text = row.names(DATA), shape = Conditions,color = col.legend.scaled)) +
            geom_point(size = cell.size, alpha = cell.transparency) +
            scale_colour_gradient(low = cell.colors[1], high = cell.colors[2], name="") +
            xlab("Dim1") +
            ylab("Dim2") +
            ggtitle(paste(MyTitle,"for (",gene,")")) +
            theme(panel.background = element_rect(fill = back.col, colour = "black"),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  legend.key = element_rect(fill = back.col))
        }
    } else {
      myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                 text = row.names(DATA), color = col.legend)) +
        geom_point(size = cell.size, alpha = cell.transparency) +
        scale_colour_gradient(low = cell.colors[1], high = cell.colors[2], name="") +
        xlab("Dim1") +
        ylab("Dim2") +
        ggtitle(paste(MyTitle,"for (",gene,")")) +
        theme_bw()
    }
  }
  # plot 3d
  if (clust.dim == 3) {
    myPLOT <- plot_ly(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3], text = row.names(DATA),
                      color = col.legend.bin, opacity = cell.transparency, marker = list(size = cell.size + 2)) %>%
      layout(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3]) %>%
      layout(plot_bgcolor = back.col) %>%
      layout(paper_bgcolor = back.col) %>%
      layout(title = MyTitle,
             scene = list(xaxis = list(title = "Dim1"),
                          yaxis = list(title = "Dim2"),
                          zaxis = list(title = "Dim3")))
    }
  }
  #######
  # plot box plot
  Yaxis1 = data.expr[[gene]]
  Yaxis = Yaxis1 + 1
  Yaxis = log2(Yaxis)
#
  if (plot.type == "boxplot") {
  myPLOT <- ggplot(data.binary, aes(x = col.legend.box, y=Yaxis)) +
    theme_bw() + theme(axis.text.x=element_text(angle=90)) +
    geom_jitter(color = box.cell.col, size = cell.size, alpha = cell.transparency) +
    ggtitle(gene) +
    geom_violin(trim=T, col = "black", alpha = cell.transparency) +
    geom_boxplot(fill = box.color,
                 col = "green",
                 notch = F,
                 outlier.shape = NA,
                 alpha = cell.transparency) +
    ylab("scaled normalized expression") +
    stat_summary(fun.y=mean, geom="point", size=2, color="blue") +
    xlab(".")
  # add p-val
  if (box.to.test == 0) {
    AvData <- x@clust.avg
    row.names(AvData) <- AvData$gene
    AvData <- AvData[,-1]
    AvData <- subset(AvData,row.names(AvData) == gene)
    box.to.test <- as.numeric(which.max(AvData))
  }
  if (box.pval == "sig.signs") {
    myPLOT <- myPLOT + stat_compare_means(label = "p.signif", ref.group = box.to.test)
  }
  if (box.pval == "sig.values") {
    myPLOT <- myPLOT + stat_compare_means(aes(label = paste0("p = ", ..p.format..)),ref.group = box.to.test)
  }
  }
  ############################ Bar plot
  mydata=cbind(data.expr,Yaxis1,col.legend.box)
  ## function to make sd
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  ### get sd data ready
  df2 <- data_summary(mydata, varname="Yaxis1",
                      groupnames=c("col.legend.box"))
  ### plot
  if (plot.type == "barplot") {
  myPLOT <- ggplot(df2, aes(x = col.legend.box, y = Yaxis1, fill = col.legend.box)) +
    stat_summary(fun.y="mean",
                 geom="bar",
                 alpha = cell.transparency,
                 show.legend = F) +
    geom_errorbar(aes(ymin=Yaxis1-sd, ymax=Yaxis1+sd), width=.2,
                  position=position_dodge(.9)) +
    ylab("avraged normalized expression") +
    xlab(".") +
    ggtitle(gene) +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
  }
  # return
  if (interactive == T) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}
