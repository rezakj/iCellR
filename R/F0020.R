#' Plot nGenes, UMIs and perecent mito
#'
#' This function takes an object of class iCellR and creates plots to see the clusters.
#' @param x An object of class iCellR.
#' @param cell.size A numeric value for the size of the cells, default = 1.
#' @param plot.type Choose between "tsne", "pca", "umap", "knetl", "diffusion", default = "tsne".
#' @param cell.color Choose cell color if col.by = "monochrome", default = "black".
#' @param back.col Choose background color, default = "black".
#' @param col.by Choose between "clusters", "conditions", "cc" (cell cycle) or "monochrome", default = "clusters".
#' @param cond.shape If TRUE the conditions will be shown in shapes.
#' @param cond.facet Show the conditions in separate plots.
#' @param anno.clust Annotate cluster names on the plot, default = TRUE.
#' @param anno.size If anno.clust is TRUE set font size, default = 3.
#' @param cell.transparency A numeric value between 0 to 1, default = 0.5.
#' @param clonotype.max Number of clonotype to plot, default = 10.
#' @param clust.dim A numeric value for plot dimensions. Choose either 2 or 3, default = 2.
#' @param interactive If TRUE an html interactive file will be made, default = TRUE.
#' @param out.name Output name for html file if interactive = TRUE, default = "plot".
#' @param angle A number to rotate the non-interactive 3D plot.
#' @param density If TRUE the density plots for PCA/tSNE second dimension will be created, default = FALSE.
#' @param static3D If TRUE a non-interactive 3D plot will be made.
#' @return An object of class iCellR.
#' @examples
#' cluster.plot(demo.obj,plot.type = "umap",interactive = FALSE)
#'
#' cluster.plot(demo.obj,plot.type = "tsne",interactive = FALSE)
#'
#' cluster.plot(demo.obj,plot.type = "pca",interactive = FALSE)
#'
#' cluster.plot(demo.obj,plot.type = "pca",col.by = "conditions",interactive = FALSE)
#'
#' cluster.plot(demo.obj,plot.type = "umap",col.by = "conditions",interactive = FALSE)
#'
#' cluster.plot(demo.obj,plot.type = "tsne",col.by = "conditions",interactive = FALSE)
#' @import RColorBrewer
#' @import scatterplot3d
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly ggplotly layout plot_ly
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @importFrom methods new
#' @importFrom stats aggregate as.dendrogram cor cor.test dist hclust p.adjust prcomp quantile sd t.test
#' @importFrom utils capture.output packageVersion read.table write.table
#' @importFrom graphics legend par plot
#' @importFrom ggplot2 ggplot theme_classic geom_segment geom_violin guide_colorbar guide_legend guides scale_color_discrete scale_colour_gradient scale_fill_gradient2 scale_x_continuous scale_y_continuous scale_y_discrete stat_summary coord_polar element_rect element_text element_blank facet_wrap scale_color_manual geom_hline geom_jitter geom_vline ylab xlab ggtitle theme_bw aes theme geom_bar geom_point geom_boxplot geom_errorbar position_dodge geom_tile geom_density geom_line
#' @export
cluster.plot <- function (x = NULL,
                          cell.size = 0.5,
                          plot.type = "tsne",
                          cell.color = "black",
                          back.col = "white",
                          col.by = "clusters",
                          cond.facet = FALSE,
                          cond.shape = FALSE,
                          anno.clust = FALSE,
                          anno.size = 4,
                          cell.transparency = 1,
                          clust.dim = 2,
                          angle = 20,
                          clonotype.max = 10,
                          density = FALSE,
                          interactive = TRUE,
                          static3D = FALSE,
                          out.name = "plot") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  if (clust.dim != 2 && clust.dim != 3) {
    stop("clust.dim should be either 2 or 3")
  }
  ###########
  # 2 dimentions
  if (clust.dim == 2) {
    if (plot.type == "tsne") {
      MyTitle = "tSNE Plot"
      DATA <- x@tsne.data
    }
    if (plot.type == "pca") {
      MyTitle = "PCA Plot"
      DATA <- x@pca.data
    }
    if (plot.type == "umap") {
      MyTitle = "UMAP Plot"
      DATA <- x@umap.data
    }
    if (plot.type == "diffusion") {
      MyTitle = "Diffusion Map Plot"
      DATA <- x@diffusion.data
    }
    if (plot.type == "knetl") {
      MyTitle = "KNetL Plot"
      DATA <- x@knetl.data
    }
  }
  # 3 dimentions
  if (clust.dim == 3) {
    if (plot.type == "tsne") {
      MyTitle = "3D tSNE Plot"
      DATA <- x@tsne.data.3d
    }
    if (plot.type == "umap") {
      MyTitle = "3D umap Plot"
      DATA <- x@umap.data
    }
    if (plot.type == "pca") {
      MyTitle = "3D PCA Plot"
      DATA <- x@pca.data
    }
    if (plot.type == "knetl") {
      MyTitle = "KNetL Plot"
      DATA <- x@knetl.data.3d
    }
    if (plot.type == "diffusion") {
      MyTitle = "3D Diffusion Map Plot"
        DATA <- x@diffusion.data
    }
  }
  # cell cycle
  if (col.by == "cc") {
    d1 <- as.character(x@stats$Phase)
    d2 <- rownames(x@stats)
    d3 <- as.data.frame(cbind(d2,d1))
    d4 <- subset(d3, d3$d2 %in% row.names(DATA))
  col.legend <- factor(d4$d1)
  }
  # conditions
  if (col.by == "conditions") {
    Cells <- rownames(DATA)
    MYConds <- as.character((unique(data.frame(do.call('rbind', strsplit(as.character(Cells),'_',fixed=TRUE)))[1]))$X1)
    if (length(MYConds) == 1) {
      stop("You need more then one condition/sample to run this")
    }
    col.legend <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
    col.legend <- factor(as.matrix(col.legend))
    }
  # clusters
  # always use hierarchical (k means changes everytime you run)
  if (col.by == "clusters") {
    if (dim(x@best.clust)[1] == 0) {
      ha <- as.data.frame(cbind(row.names(DATA),rep(0,dim(DATA)[1])))
      colnames(ha) <- c("rows","clusters")
      row.names(ha) <- ha$rows
      ha <- as.data.frame(as.matrix(ha)[,-1])
      colnames(ha) <- "clusters"
#      ha$clusters <- as.numeric(ha$clusters)
      x@best.clust <- ha
      col.legend <- factor(x@best.clust$clusters)
    } else {
      col.legend <- factor(x@best.clust$clusters)
    }
  }
  ## clonotype
  if (col.by == "clonotype") {
    if (is.null(x@vdj.data)) {
      stop("Clonotype data is missing use add.vdj function to load them.")
    } else {
      colono <- unique(x@vdj.data[1:2])
      row.names(colono) <- cell.barcodes
      colono$raw_clonotype_id <- gsub("clonotype"," ", colono$raw_clonotype_id)
      colono <- colono[1]
      colnames(colono) <- c("Clonotypes")
      colonoData <- merge(DATA,colono, by="row.names", all.x=TRUE, all.y=FALSE)
      colonoData$Clonotypes[as.numeric(colonoData$Clonotypes) > clonotype.max] <- "remaining"
      colonoData$Clonotypes[is.na(colonoData$Clonotypes)] <- "not.determined"
      colonoData$Clonotypes <- gsub( " ", "", colonoData$Clonotypes)
      col.legend <- factor(colonoData$Clonotypes)
      X <- (unique(colonoData$Clonotypes))
      X1 <- as.character(sort(suppressWarnings(as.numeric(X))))
      X2 <- c(X1,"remaining","not.determined")
      col.legend <- factor(col.legend, levels = X2)
    }
  }
  # monochrome
  if (col.by == "monochrome") {
    if (clust.dim == 2) {
    myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                               text = row.names(DATA))) +
      geom_point(size = cell.size, col = cell.color, alpha = cell.transparency) +
      guides(colour = guide_legend(override.aes = list(size=5))) +
      xlab("Dim1") +
      ylab("Dim2") +
      ggtitle(MyTitle) +
      scale_color_discrete(name="") +
      theme(panel.background = element_rect(fill = back.col, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = back.col)) + theme_bw()
    }
    if (clust.dim == 3) {
      myPLOT <- plot_ly(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3], text = row.names(DATA),
              opacity = cell.transparency, marker = list(size = cell.size + 2, color = cell.color)) %>%
        layout(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3]) %>%
        layout(plot_bgcolor = back.col) %>%
        layout(paper_bgcolor = back.col) %>%
        layout(title = MyTitle,
             scene = list(xaxis = list(title = "Dim1"),
                          yaxis = list(title = "Dim2"),
                          zaxis = list(title = "Dim3")))
      }
    }
# plot 2d
  if (clust.dim == 2) {
    if (cond.shape == FALSE) {
      if (static3D == FALSE) {
        myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                   text = row.names(DATA), color = col.legend)) +
          geom_point(size = cell.size, alpha = cell.transparency) +
          guides(colour = guide_legend(override.aes = list(size=5))) +
          xlab("Dim1") +
          ylab("Dim2") +
          ggtitle(MyTitle) +
          scale_color_discrete(name="") +
          theme(panel.background = element_rect(fill = back.col, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.key = element_rect(fill = back.col))
############# Annotation
        if (col.by == "clusters") {
          if (anno.clust == TRUE) {
            df <- cbind(x@best.clust, DATA)
            cords <- aggregate(df[, 2:3], list(df$clusters), mean)
            MYX=cords$V1
            MYY=cords$V2
            MYZ=cords$Group.1
            myPLOT <- myPLOT +  annotate("text",
                                         x = MYX,
                                         y = MYY,
                                         label = MYZ,
                                         colour = "black",
                                         size = anno.size)
          }
        }
#############
      }
      if (interactive == TRUE) {
        myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                   text = row.names(DATA), color = col.legend)) +
          geom_point(size = cell.size, alpha = cell.transparency) +
          guides(colour = guide_legend(override.aes = list(size=5))) +
          xlab("Dim1") +
          ylab("Dim2") +
          ggtitle(MyTitle) +
          scale_color_discrete(name="") +
          theme_bw()
      }
    }
      if (cond.shape == TRUE) {
        conds.sh <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
        cond.shape <- factor(as.matrix(conds.sh))
        if (interactive == FALSE) {
          myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                     text = row.names(DATA), color = col.legend, shape = cond.shape)) +
            geom_point(size = cell.size, alpha = cell.transparency) +
            guides(colour = guide_legend(override.aes = list(size=5))) +
            xlab("Dim1") +
            ylab("Dim2") +
            ggtitle(MyTitle) +
            scale_color_discrete(name="") +
            theme(panel.background = element_rect(fill = back.col, colour = "black"),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  legend.key = element_rect(fill = back.col))
        }
        if (interactive == TRUE) {
          myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                     text = row.names(DATA), color = col.legend, shape = cond.shape)) +
            geom_point(size = cell.size, alpha = cell.transparency) +
            guides(colour = guide_legend(override.aes = list(size=5))) +
            xlab("Dim1") +
            ylab("Dim2") +
            ggtitle(MyTitle) +
            scale_color_discrete(name="") +
            theme_bw()
      }
      }
    if (cond.facet == TRUE) {
      if (interactive == FALSE) {
      conds.sh <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
      cond.shape <- factor(as.matrix(conds.sh))
      myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                 text = row.names(DATA), color = col.legend)) +
        geom_point(size = cell.size, alpha = cell.transparency) +
        guides(colour = guide_legend(override.aes = list(size=5))) +
        xlab("Dim1") +
        ylab("Dim2") +
        ggtitle(MyTitle) +
        scale_color_discrete(name="") +
        theme(panel.background = element_rect(fill = back.col, colour = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.key = element_rect(fill = back.col)) + facet_wrap(cond.shape)
      }
    }
  }
# plot 3d
  if (clust.dim == 3) {
    if (static3D == FALSE) {
#      DATAann <- as.data.frame(x@cluster.data$Best.partition)
      DATAann <- as.data.frame(x@best.clust)
      A = (row.names(DATAann))
      colnames(DATAann) <- ("clusters")
      B = as.character(as.matrix(DATAann[1]))
      ANNOT <- paste(A,"clust",B,sep="_")
    myPLOT <- plot_ly(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3], text = ANNOT,
                      color = col.legend, opacity = cell.transparency, marker = list(size = cell.size + 2)) %>%
      layout(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3]) %>%
      layout(plot_bgcolor = back.col) %>%
      layout(paper_bgcolor = back.col) %>%
      layout(title = MyTitle,
           scene = list(xaxis = list(title = "Dim1"),
                        yaxis = list(title = "Dim2"),
                        zaxis = list(title = "Dim3")))
  } else {
    # colors
#    col.legend <- factor(x@best.clust$clusters)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    colors <- c(colors[1:3],colors[5:30])
    if (col.by == "clusters") {
      cols <- colors[as.numeric(x@best.clust$clusters)]
      col.legend <- as.factor(x@best.clust$clusters)

    }
    if (col.by == "conditions") {
      col.legend <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
      col.legend <- factor(as.matrix(col.legend))
      cols <- colors[col.legend]
    }
    #
    scatterplot3d (x = DATA[,2], y = DATA[,3], z = DATA[,1],
                color = cols,
                pch = 19,
                xlab = "Dim2", ylab = "Dim3",zlab = "Dim1",
                main = MyTitle,
                grid = TRUE,
                box = TRUE,
                scale.y = 1,
                angle = angle,
                mar = c(3,3,3,6)+0.1,
                cex.axis = 0.6,
                cex.symbols = cell.size,
                highlight.3d = FALSE)
    # legend
    legend("topright", legend = levels(col.legend),
           col =  colors,
           pch = 19,
           inset = -0.1,
           xpd = TRUE,
           horiz = FALSE)
#    scatter3D(x = DATA[,2], y = DATA[,3], z = DATA[,1],
#              colvar = NULL,
#              col = col.legend,
#              pch = 19,
#              cex = 1,
#              colkey = T,
#              col.panel ="steelblue",
#              col.grid = "darkblue",
#              ticktype = "detailed",
#              bty ="g",
#              d = 2,
#             phi = 10,
#              theta = 15)
    myPLOT = "plot made"
  }
}
  # density plot
  if (density == TRUE) {
    myPLOT <- ggplot(DATA, aes(DATA[,2], fill = col.legend)) +
      geom_density(alpha=cell.transparency) +
      xlab("Dim2") +
      scale_color_discrete(name="") + theme_bw()
  }
# return
  if (interactive == TRUE) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
    }
}

