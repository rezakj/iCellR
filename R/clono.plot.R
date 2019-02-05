#' Make scatter, box and bar plots for genes
#'
#' This function takes an object of class scSeqR and provides plots for genes.
#' @param x An object of class scSeqR.
#' @param gene A gene name to be plotted.
#' @param box.to.test A cluster number so that all the boxes in the box plot would be compared to. If set to "0" the cluster with the highest avrage would be choosen, defult = 0.
#' @param box.pval Choose from "sig.values" and "sig.signs". If set to "sig.signs" p values would be replaced with signs ("na", "*", "**", "***"), defult = "sig.signs".
#' @param plot.data.type Choose from "tsne" and "pca", defult = "tsne".
#' @param clust.dim 2 for 2D plots and 3 for 3D plots, defult = 2.
#' @param col.by Choose from "clusters" and "conditions", defult = "clusters".
#' @param plot.type Choose from "scatterplot", "boxplot" and "barplot", defult = "scatterplot".
#' @param cell.size A number for the size of the points in the plot, defult = 1.
#' @param cell.colors Colors for heat mapping the points in "scatterplot", defult = c("gray","red").
#' @param box.cell.col A color for the points in the box plot, defult = "black".
#' @param box.color A color for the boxes in the "boxplot", defult = "red".
#' @param box.line.col A color for the lines around the "boxplot", defult = "green".
#' @param back.col A color for the plot background, defult = "black".
#' @param cell.transparency Color transparency for points in "scatterplot" and "boxplot", defult = 0.5.
#' @param interactive If set to TRUE an intractive HTML file will be created, defult = TRUE.
#' @param out.name If "interactive" is set to TRUE, the out put name for HTML, defult = "plot".
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' clono.plot(my.obj,
#'             plot.data.type = "tsne",
#'             clono = 1,
#'             clust.dim = 2,
#'             cell.size = 1,
#'             cell.colors = c("red","gray"),
#'             cbox.cell.col = "black",
#'             back.col = "white",
#'             interactive = T,
#'             cell.transparency = 0.5,
#'             out.name = "tSNE_3D_clusters")
#' }
#' @import ggpubr
#' @export
clono.plot <- function (x = NULL,
                       plot.data.type = "tsne",
                       clono = 1,
                       clust.dim = 2,
                       cell.size = 1,
                       cell.colors = c("red","gray"),
                       box.cell.col = "black",
                       back.col = "white",
                       cell.transparency = 0.5,
                       interactive = TRUE,
                       out.name = "plot") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
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
  }
  ## clonotype
      colono <- unique(x@vdj.data[1:2])
      row.names(colono) <- colono$barcode
      colono$raw_clonotype_id <- gsub("clonotype"," ", colono$raw_clonotype_id)
      colono <- colono[1]
      colnames(colono) <- c("Clonotypes")
      colonoData <- merge(DATA,colono, by="row.names", all.x=T, all.y=F)
      colonoData$Clonotypes <- gsub( " ", "", colonoData$Clonotypes)
      colonoData$Clonotypes[is.na(colonoData$Clonotypes)] <- "NA"
      colonoData$Clonotypes[colonoData$Clonotypes != clono] <- "NA"
      colonoData$Clonotypes[colonoData$Clonotypes == clono] <- clono
      DATA <- colonoData
      row.names(DATA) <- DATA$Row.names
      DATA <- DATA[,-1]
      DATA <- (DATA[order(DATA$Clonotypes, decreasing = T),])
      clonotype <- factor(DATA$Clonotypes)
    # plot 2d
    if (clust.dim == 2) {
      if (interactive == F) {
        myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],col=clonotype,
                                   text = row.names(DATA))) +
          geom_point(size = cell.size, alpha = cell.transparency) +
          xlab("Dim1") +
          ylab("Dim2") +
          ggtitle(paste(MyTitle)) +
          scale_color_manual(values = cell.colors) +
          theme(panel.background = element_rect(fill = back.col, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.key = element_rect(fill = back.col))
      } else {
        myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                   text = row.names(DATA), color = clonotype)) +
          geom_point(size = cell.size, alpha = cell.transparency) +
          scale_color_manual(values = cell.colors) +
          xlab("Dim1") +
          ylab("Dim2") +
          ggtitle(MyTitle) +
          theme_bw()
      }
    }
    # plot 3d
    if (clust.dim == 3) {
      myPLOT <- plot_ly(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3], text = row.names(DATA),
                        color = clonotype, opacity = cell.transparency, marker = list(size = cell.size + 2)) %>%
        layout(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3]) %>%
        layout(plot_bgcolor = back.col) %>%
        layout(paper_bgcolor = back.col) %>%
        layout(title = MyTitle,
               scene = list(xaxis = list(title = "Dim1"),
                            yaxis = list(title = "Dim2"),
                            zaxis = list(title = "Dim3")))
    }
  # return
  if (interactive == T) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}
