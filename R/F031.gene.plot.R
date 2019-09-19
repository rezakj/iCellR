#' Make scatter, box and bar plots for genes
#'
#' This function takes an object of class iCellR and provides plots for genes.
#' @param x An object of class iCellR.
#' @param gene A gene name to be plotted.
#' @param box.to.test A cluster number so that all the boxes in the box plot would be compared to. If set to "0" the cluster with the highest avrage would be choosen, default = 0.
#' @param box.pval Choose from "sig.values" and "sig.signs". If set to "sig.signs" p values would be replaced with signs ("na", "*", "**", "***"), default = "sig.signs".
#' @param plot.data.type Choose between "tsne", "pca", "umap", "diffusion", "pseudo.A" and "pseudo.B", default = "tsne".
#' @param clust.dim 2 for 2D plots and 3 for 3D plots, default = 2.
#' @param col.by Choose from "clusters" and "conditions", default = "clusters".
#' @param data.type Choose from "main" or "imputed", default = "main".
#' @param scaleValue Scale the colors, default = FALSE.
#' @param min.scale If scaleValue = TRUE, set a number for min, default = -2.5.
#' @param max.scale If scaleValue = TRUE, set a number for max, default = 2.5.
#' @param cond.shape If TRUE the conditions will be shown in shapes.
#' @param conds.to.plot Choose the conditions you want to see in the plot, default = NULL (all conditions).
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
#' gene.plot(demo.obj, gene = "CD74",interactive = FALSE)
#'
#' gene.plot(demo.obj, gene = "CD74",plot.data.type = "umap",interactive = FALSE)
#'
#' gene.plot(demo.obj, gene = "CD74",
#'           plot.data.type = "umap",
#'           interactive = FALSE,
#'           plot.type = "barplot")
#'
#' gene.plot(demo.obj, gene = "CD74",
#'           plot.data.type = "umap",
#'           interactive = FALSE,
#'           plot.type = "boxplot")
#'
#' @importFrom ggpubr stat_compare_means
#' @import plyr
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly ggplotly layout plot_ly
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @importFrom methods new
#' @importFrom stats aggregate as.dendrogram cor cor.test dist hclust p.adjust prcomp quantile sd t.test
#' @importFrom utils capture.output packageVersion read.table write.table
#' @importFrom graphics legend par plot
#' @importFrom ggplot2 ggplot geom_segment geom_violin guide_colorbar guide_legend guides scale_color_discrete scale_colour_gradient scale_fill_gradient2 scale_x_continuous scale_y_continuous scale_y_discrete stat_summary coord_polar element_rect element_text element_blank facet_wrap scale_color_manual geom_hline geom_jitter geom_vline ylab xlab ggtitle theme_bw aes theme geom_bar geom_point geom_boxplot geom_errorbar position_dodge geom_tile geom_density geom_line
#' @export
gene.plot <- function (x = NULL,
                       gene = "NULL",
                       cond.shape = FALSE,
                       conds.to.plot = NULL,
                       data.type = "main",
                       box.to.test = 0,
                       box.pval = "sig.signs",
                       plot.data.type = "tsne",
                       scaleValue = FALSE,
                       min.scale = -2.5,
                       max.scale = 2.5,
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
    if (plot.type == "pseudo.A") {
      MyTitle = "Pseudo Map A Plot"
      DATA <- x@pseudo.mapA
    }
    if (plot.type == "pseudo.B") {
      MyTitle = "Pseudo Map B Plot"
      DATA <- x@pseudo.mapB
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
    MyConds = col.legend.box
  }
  # clusters
  if (col.by == "clusters") {
    if (is.null(x@best.clust)) {
      stop("Clusters are not assigend yet")
    } else {
      col.legend.box <- x@best.clust
      col.legend.box <- factor(x@best.clust$clusters)
#      col.legend.box$clusters <- sub("^", "cl.",col.legend.box$clusters)
#      col.legend.box <- factor(col.legend.box$clusters)
    }
  }
  # get the gene from the main data
  sub.data <- subset(DATAmain,rownames(DATAmain) == gene)
  data.t <- t(sub.data)
  data.expr <- as.data.frame(data.t)
###### make binary
#  data.binary <- data.t > 0
#  data.binary <- as.data.frame(data.binary)
#  col.legend.bin = data.binary[[gene]]
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
# fix color for ADTs
  if (scaleValue == FALSE) {
      if ( length(grep("^ADT_", gene, value = TRUE)) == 1) {
        Lo3=(quantile(col.legend,0.05)) # 1
        Lo2=(quantile(col.legend,0.10)) # 2
        Lo1=(quantile(col.legend,0.25)) # 3
        MID=(quantile(col.legend,0.50)) # 4
        Up1=(quantile(col.legend,0.75)) # 5
        Up2=(quantile(col.legend,0.90)) # 6
        Up3=(quantile(col.legend,0.95)) # 7
        col.legend <- replace(col.legend, col.legend > Up3, Up3)
        col.legend <- replace(col.legend, col.legend > Up2 & col.legend < Up3 ,Up2)
        col.legend <- replace(col.legend, col.legend > Up1 & col.legend < Up2 ,Up1)
        col.legend <- replace(col.legend, col.legend > MID & col.legend < Up1 ,MID)
        col.legend <- replace(col.legend, col.legend > Lo1 & col.legend < MID ,Lo1)
        col.legend <- replace(col.legend, col.legend < Lo2 ,Lo1)
      }
  }
#      if ( length(grep("^ADT_", gene, value = T)) == 0) {
#        col.legend = col.legend
#      }
###
  if (scaleValue == TRUE) {
    col.legend <- scale(col.legend)
    col.legend <- FixScale(mydata = col.legend, min = min.scale, max = max.scale)
  }
###
# fix the dataframe
MyConds <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
MyConds <- factor(as.matrix(MyConds))
DATA <- cbind(DATA, Expression = col.legend, Clusters = col.legend.box, Conditions = MyConds)
if (!is.null(conds.to.plot)) {
DATA <- subset(DATA, DATA$Conditions %in% conds.to.plot)
}
###
  if (plot.type == "scatterplot") {
  # plot 2d
#    Conditions = col.legend.box
  if (clust.dim == 2) {
      if (cond.shape == FALSE) {
        myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                   text = row.names(DATA), color = Expression)) +
          geom_point(size = cell.size, alpha = cell.transparency) +
          scale_colour_gradient(low = cell.colors[1], high = cell.colors[2], name="") +
          xlab("Dim1") +
          ylab("Dim2") +
          ggtitle(paste(MyTitle,"for (",gene,")")) +
          theme(panel.background = element_rect(fill = back.col, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.key = element_rect(fill = back.col))
      }
        if (cond.shape == TRUE) {
          myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                     text = row.names(DATA), shape = Conditions,color = Expression)) +
            geom_point(size = cell.size, alpha = cell.transparency) +
            scale_colour_gradient(low = cell.colors[1], high = cell.colors[2], name="") +
            xlab("Dim1") +
            ylab("Dim2") +
            ggtitle(paste(MyTitle,"for (",gene,")")) +
            theme(panel.background = element_rect(fill = back.col, colour = "black"),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  legend.key = element_rect(fill = back.col))
    } #else {
#      myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
#                                 text = row.names(DATA), color = Expression)) +
#        geom_point(size = cell.size, alpha = cell.transparency) +
#        scale_colour_gradient(low = cell.colors[1], high = cell.colors[2], name="") +
#        xlab("Dim1") +
#        ylab("Dim2") +
#        ggtitle(paste(MyTitle,"for (",gene,")")) +
#        theme_bw()
#    }
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
#  Yaxis1 = data.expr[[gene]]
#  Yaxis = Yaxis1 + 1
#  Yaxis = log2(Yaxis)
#
  if (plot.type == "boxplot") {
    if (cond.shape == FALSE) {
  myPLOT <- ggplot(DATA, aes(x = Clusters, y=log2(Expression + 1))) +
    theme_bw() + theme(axis.text.x=element_text(angle=90)) +
    geom_jitter(color = box.cell.col, size = cell.size, alpha = cell.transparency) +
    ggtitle(gene) +
    geom_violin(trim=TRUE, col = "black", alpha = cell.transparency) +
    geom_boxplot(fill = box.color,
                 col = "green",
                 notch = FALSE,
                 outlier.shape = NA,
                 alpha = cell.transparency) +
    ylab("scaled normalized expression") +
    stat_summary(fun.y=mean, geom="point", size=2, color="blue") +
    xlab(".")
    }
    if (cond.shape == TRUE) {
#      Conditions = DATA$Conditions
      myPLOT <- ggplot(DATA, aes(x = Clusters, y=log2(Expression + 1))) +
        theme_bw() + theme(axis.text.x=element_text(angle=90)) +
        geom_jitter(color = box.cell.col, size = cell.size, alpha = cell.transparency) +
        ggtitle(gene) +
        geom_violin(trim=TRUE, col = "black", alpha = cell.transparency) +
        geom_boxplot(fill = box.color,
                     col = "green",
                     notch = FALSE,
                     outlier.shape = NA,
                     alpha = cell.transparency) +
        ylab("scaled normalized expression") +
        stat_summary(fun.y=mean, geom="point", size=2, color="blue") +
        xlab(".")
      myPLOT <- myPLOT + facet_wrap(~ Conditions)
    }
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
#  mydata=cbind(data.expr,Yaxis1,col.legend.box)
  ## function to make sd
  data_summary <- function(data, varname, groupnames){
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
  if (plot.type == "barplot") {
  if (cond.shape == FALSE) {
    df2 <- data_summary(DATA, varname="Expression",
                        groupnames=c("Clusters"))
    myPLOT <- ggplot(df2, aes(x = Clusters, y = Expression, fill = Clusters)) +
      geom_errorbar(aes(ymin=Expression, ymax=Expression+sd), width=.2,
                    position=position_dodge(.9)) +
      stat_summary(fun.y="mean",
                   geom="bar",
                   alpha = cell.transparency,
                   show.legend = FALSE) +
      ylab("avraged normalized expression") +
      xlab(".") +
      ggtitle(gene) +
      theme_bw() + theme(axis.text.x=element_text(angle=90))
  }
  if (cond.shape == TRUE) {
    df2 <- data_summary(DATA, varname="Expression",
                        groupnames=c("Conditions","Clusters"))
    myPLOT <- ggplot(df2, aes(x = Clusters, y = Expression, fill = Clusters)) +
      geom_errorbar(aes(ymin=Expression, ymax=Expression+sd), width=.2,
                    position=position_dodge(.9)) +
      stat_summary(fun.y="mean",
                   geom="bar",
                   alpha = cell.transparency,
                   show.legend = FALSE) +
      ylab("avraged normalized expression") +
      xlab(".") +
      ggtitle(gene) +
      theme_bw() + theme(axis.text.x=element_text(angle=90)) + facet_wrap(~ Conditions)
   }
  }
  # return
  if (interactive == TRUE) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}
