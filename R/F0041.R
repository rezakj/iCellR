#' Make 2D and 3D scatter plots for clonotypes.
#'
#' This function takes an object of class iCellR and provides plots for clonotypes.
#' @param x An object of class iCellR.
#' @param clono A clonotype name to be plotted, default = NULL.
#' @param clonotype.column The column which has the clonotype IDs, default =  2.
#' @param barcode.column The column which has the barcode IDs, default = 1.
#' @param plot.data.type Choose from "tsne" and "pca", default = "tsne".
#' @param conds.to.plot Choose one condition you want to see in the plot, default = NULL (all conditions).
#' @param clust.dim 2 for 2D plots and 3 for 3D plots, default =  2.
#' @param box.cell.col Choose a color for box default =  "black".
#' @param cell.size A number for the size of the points in the plot, default = 1.
#' @param cell.colors Colors for heat mapping the points in "scatterplot", default = c("gray","red").
#' @param back.col A color for the plot background, default = "black".
#' @param cell.transparency Color transparency for points, default = 0.5.
#' @param interactive If set to TRUE an intractive HTML file will be created, default = TRUE.
#' @param out.name If "interactive" is set to TRUE, the out put name for HTML, default = "plot".
#' @return An object of class iCellR.
#' @export
clono.plot <- function (x = NULL,
                       plot.data.type = "tsne",
                       clonotype.column = 1,
                       barcode.column = 2,
                       clono = NULL,
                       conds.to.plot = NULL,
                       clust.dim = 2,
                       cell.size = 1,
                       cell.colors = c("red","gray"),
                       box.cell.col = "black",
                       back.col = "white",
                       cell.transparency = 1,
                       interactive = TRUE,
                       out.name = "plot") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ###
  if (is.null(clono)) {
    stop("Please provide a clonotype id (example: S1_clonotype1)")
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
    if (plot.data.type == "knetl") {
      MyTitle = "KNetL Plot"
      DATA <- x@knetl.data
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
    if (plot.data.type == "knetl") {
      MyTitle = "KNetL Plot"
      DATA <- x@knetl.data.3d
    }
    if (plot.data.type == "diffusion") {
      MyTitle = "Diffusion Map Plot"
      DATA <- x@diffusion.data
    }
  }
  #######
  MyRows <- rownames(DATA)
  rownames(DATA) <- gsub("-",".",MyRows)
  ## clonotype

  Bardata <- as.character(as.matrix(x@vdj.data[barcode.column]))
  Colodata <- as.character(as.matrix(x@vdj.data[clonotype.column]))
      colono <- as.data.frame(cbind(Bardata,Colodata))
      ######
      colnames(colono) <- c("barcode","Clonotypes")
      colono <- unique(colono)
      colono$barcode <- gsub("-",".",colono$barcode)
      colono <- subset(colono, colono$Clonotypes == clono)
#row.names(colono) <- as.character(colono$barcode)
      ####### Conditions
           if (!is.null(conds.to.plot)) {
             col.legend <- data.frame(do.call('rbind', strsplit(as.character(colono$barcode),'_',fixed=TRUE)))[1]
             colnames(col.legend) <- "Conditions"
             colono$Conditions <- col.legend
             colono <- subset(colono, colono$Conditions == conds.to.plot)
           }
#     head(colono)
#      colono <- unique(x@vdj.data[1:2])
#      cell.barcodes <- colono$barcode
#      cell.barcodes <- gsub("-",".",cell.barcodes)
#      row.names(colono) <- cell.barcodes

#           colono$raw_clonotype_id <- gsub("clonotype"," ", colono$raw_clonotype_id)
#      colono <- colono[1]

#       if (dim(do.call('rbind', strsplit(as.character(rownames(colono)),'_',fixed=TRUE)))[2] == 1) {
#        if (dim(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[2] == 2) {
#          rownames(DATA) <- as.character(as.matrix(as.data.frame((do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE))))[2]))
#        }
#       }
############## merge
#     clono = "S1_clonotype1"
#############
     DATA1 <- subset(DATA, rownames(DATA) %in% colono$barcode)
     DATA1$MyCol <- cell.colors[1]
     DATA2 <- subset(DATA, !rownames(DATA) %in% colono$barcode)
     DATA2$MyCol <- cell.colors[2]
     DATA <- rbind(DATA2,DATA1)
#############
     MyCol <- factor(DATA$MyCol)
    # plot 2d
    if (clust.dim == 2) {
      if (interactive == FALSE) {
        myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                   text = row.names(DATA))) +
          geom_point(size = cell.size, alpha = cell.transparency, col=MyCol) +
          xlab("Dim1") +
          ylab("Dim2") +
          ggtitle(paste(MyTitle, "(clonotype",clono, ")")) +
          theme(legend.position = "none") +
          scale_color_manual(values = cell.colors) +
          theme(panel.background = element_rect(fill = back.col, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.key = element_rect(fill = back.col))
      } else {
        myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                   text = row.names(DATA))) +
          geom_point(size = cell.size, alpha = cell.transparency, col=MyCol) +
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
                        color = MyCol, opacity = cell.transparency, marker = list(size = cell.size + 2)) %>%
        layout(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3]) %>%
        layout(plot_bgcolor = back.col) %>%
        layout(paper_bgcolor = back.col) %>%
        layout(title = MyTitle,
               scene = list(xaxis = list(title = "Dim1"),
                            yaxis = list(title = "Dim2"),
                            zaxis = list(title = "Dim3")))
    }
  # return
  if (interactive == TRUE) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}
