#' Plot nGenes, UMIs and perecent mito, genes, clusters and more on spatial image
#'
#' This function takes an object of class iCellR and creates spatial plots.
#' @param x An object of class iCellR.
#' @param cell.size A numeric value for the size of the cells, default = 1.
#' @param cell.colors Colors for heat mapping the points in "scatterplot", default = c("gray","red").
#' @param col.by Choose between "clusters", "mt","UMIs","nGenes", "cc" (cell cycle) or "gene", default = "clusters".
#' @param gene Gene name/names to be plotted, if col.by = "gene".
#' @param conds.to.plot Choose the conditions you want to see in the plot, default = NULL (all conditions).
#' @param data.type Choose from "main" or "imputed", default = "main".
#' @param back.col A color for the plot background, default = "black".
#' @param scaleValue Scale the colors, default = FALSE.
#' @param min.scale If scaleValue = TRUE, set a number for min, default = -2.5.
#' @param max.scale If scaleValue = TRUE, set a number for max, default = 2.5.
#' @param anno.clust Annotate cluster names on the plot, default = TRUE.
#' @param anno.size If anno.clust is TRUE set font size, default = 3.
#' @param anno.col If anno.clust is TRUE set color, default = "white".
#' @param cell.transparency Color transparency for points in "scatterplot" and "boxplot", default = 1.
#' @param interactive If set to TRUE an interactive HTML file will be created, default = TRUE.
#' @param out.name If "interactive" is set to TRUE, the out put name for HTML, default = "plot".
#' @return An object of class iCellR.
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
spatial.plot <- function (x = NULL,
                          cell.size = 1,
                          cell.colors = c("gray","red"),
                          back.col = "black",
                          col.by = "clusters",
                          conds.to.plot = NULL,
                          gene = NULL,
                          data.type = "main",
                          scaleValue = TRUE,
                          min.scale = 0,
                          max.scale = 2.5,
                          anno.clust = FALSE,
                          anno.size = 4,
                          anno.col = "white",
                          cell.transparency = 1,
                          interactive = TRUE,
                          out.name = "plot") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ###########
     DATA <- x@spatial.data
     DATA$V5 <- (DATA$V5)*-1
########
      Cells <- row.names(DATA)
      MYConds <- as.character((unique(data.frame(do.call('rbind', strsplit(as.character(Cells),'_',fixed=TRUE)))[1]))$X1)
      if (length(MYConds) == 1) {
        MYconds <- rep(NA,dim(DATA)[1])
        DATA <- cbind(DATA,MYconds)
      }
      if(length(MYConds) > 1) {
        MYconds <- as.character(as.matrix(as.data.frame(
          do.call('rbind', strsplit(as.character(Cells),'_',fixed=TRUE)))[1]))
        DATA <- cbind(DATA,MYconds)
      }
   # clusters
  # always use hierarchical (k means changes everytime you run)
  if (col.by == "clusters") {
      MyTitle="clusters"
      data <- merge(x@best.clust,DATA,by="row.names", all.x=F, all.y=F)
      data$clusters <- as.factor(data$clusters)
      }
if (col.by == "cc") {
        MyTitle="Cell Cycle"
        NewDat <- as.data.frame(cbind(row.names(x@stats), x@stats$Phase))
        row.names(NewDat) <- NewDat$V1
        NewDat <- NewDat[2]
        colnames(NewDat) <- "clusters"
        data <- merge(NewDat,DATA,by="row.names", all.x=F, all.y=F)
        }
      if (col.by == "nGenes") {
        MyTitle="Number of genes"
        NewDat <- as.data.frame(cbind(row.names(x@stats), x@stats$nGenes))
        row.names(NewDat) <- NewDat$V1
        NewDat <- NewDat[2]
        colnames(NewDat) <- "clusters"
        data <- merge(NewDat,DATA,by="row.names", all.x=F, all.y=F)
        data$clusters <- as.numeric(data$clusters)
      }
      if (col.by == "UMIs") {
        MyTitle="Number of UMIs"
        NewDat <- as.data.frame(cbind(row.names(x@stats), x@stats$UMIs))
        row.names(NewDat) <- NewDat$V1
        NewDat <- NewDat[2]
        colnames(NewDat) <- "clusters"
        data <- merge(NewDat,DATA,by="row.names", all.x=F, all.y=F)
        data$clusters <- as.numeric(data$clusters)
      }
      if (col.by == "mt") {
        MyTitle="mito.percent"
        NewDat <- as.data.frame(cbind(row.names(x@stats), x@stats$mito.percent))
        row.names(NewDat) <- NewDat$V1
        NewDat <- NewDat[2]
        colnames(NewDat) <- "clusters"
        data <- merge(NewDat,DATA,by="row.names", all.x=F, all.y=F)
        data$clusters <- as.numeric(data$clusters)
      }
      if (col.by == "gene") {
        if (data.type == "main") {
          DATAmain <- x@main.data
        }
        if (data.type == "imputed") {
          DATAmain <- x@imputed.data
        }
        if (is.null(gene)) {
          stop("There is no gene name provided. Please provide gene/genes to plot")
        }
        AllGenes = row.names(DATAmain)
        absent = which((gene %in% AllGenes) == FALSE)
        absentgenes = gene[absent]
        gene <- setdiff(gene,absentgenes)
        #####
        if(length(absentgenes) != 0)
        {
          absentgenes = paste(absentgenes, collapse=",")
          ToPrint <- paste("WARNING:",absentgenes, "not available in your data", sep=" ")
          message(ToPrint)
          presentgenes = paste(gene, collapse=",")
          ToPrint2 <- paste("WARNING:","plotting",presentgenes, sep=" ")
          message(ToPrint2)
        }
        # get the gene from the main data
        sub.data <- subset(DATAmain,rownames(DATAmain) %in% gene)
        data.t <- t(sub.data)
        MyTitle=gene
        if(1 < dim(data.t)[2]) {
          data.t <- rowSums(data.t)
          MyTitle="Multiple genes"
        }
        data.expr <- as.data.frame(data.t)
        colnames(data.expr) <- "clusters"
        data <- merge(data.expr,DATA,by="row.names", all.x=F, all.y=F)
        data$clusters <- as.numeric(data$clusters)
      }
############ scaling
      # fix scale
      FixScale <- function (mydata, min, max){
        Mydat <- mydata
        Mydat[Mydat > max] <- max
        Mydat[Mydat < min] <- min
        return(Mydat)
      }
      ###
      if (scaleValue == TRUE) {
        if (is.numeric(data$clusters)) {
        col.legend = scale(data$clusters)
        data$clusters <- FixScale(mydata = col.legend, min = min.scale, max = max.scale)
        }
      }
      ###
####### Conditions
      if (!is.null(conds.to.plot)) {
        data <- subset(data, data$MYconds %in% conds.to.plot)
      }
################# Fix V5
#     data$V5 <- data$V5 + ((min(data$V5)) * -1)
############################# plot
      if (!is.numeric(data$clusters)) {
      if(length(MYConds) == 1) {
        myPLOT <- ggplot(data, aes(V6, y = V5,
                                   text = row.names(data), color = clusters)) +
          geom_point(size = cell.size, alpha = cell.transparency) +
          guides(colour = guide_legend(override.aes = list(size=5))) +
          scale_color_discrete(name="") +
          ggtitle(MyTitle) +
          theme(panel.background = element_rect(fill = back.col, colour = back.col),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.key = element_rect(fill = back.col)) +
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) +
          theme(axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank())
        ############# Annotation
          if (anno.clust == TRUE) {
            cords <- aggregate(data[, 6:7], list(data$clusters), mean)
#            MYX=(cords$V5) * -1
#            MYY=(cords$V6) * -1
            MYX=(cords$V6)
            MYY=(cords$V5)
            MYZ=cords$Group.1
            myPLOT <- myPLOT +  annotate("text",
                                         x = MYX,
                                         y = MYY,
                                         label = MYZ,
                                         colour = anno.col,
                                         size = anno.size)
          }
      }
      if(length(MYConds) > 1) {
          if(length(conds.to.plot) == 1) {
            myPLOT <- ggplot(data, aes(V6, y = V5,
                                       text = row.names(data), color = clusters)) +
              geom_point(size = cell.size, alpha = cell.transparency) +
              guides(colour = guide_legend(override.aes = list(size=5))) +
              scale_color_discrete(name="") +
              ggtitle(MyTitle) +
              theme(panel.background = element_rect(fill = back.col, colour = back.col),
                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    legend.key = element_rect(fill = back.col)) +
              theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()) +
              theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
            ############# Annotation
              if (anno.clust == TRUE) {
                cords <- aggregate(data[, 6:7], list(data$clusters), mean)
                MYX=(cords$V6)
                MYY=(cords$V5)
                MYZ=cords$Group.1
                myPLOT <- myPLOT +  annotate("text",
                                             x = MYX,
                                             y = MYY,
                                             label = MYZ,
                                             colour = anno.col,
                                             size = anno.size)
            }
            #######
          }
        if(length(conds.to.plot) > 1) {
          if(length(conds.to.plot) > 1) {
            myPLOT <- ggplot(data, aes(V6, y = V5,
                                       text = row.names(data), color = clusters)) +
              geom_point(size = cell.size, alpha = cell.transparency) +
              guides(colour = guide_legend(override.aes = list(size=5))) +
              scale_color_discrete(name="") +
              ggtitle(MyTitle) +
              theme(panel.background = element_rect(fill = back.col, colour = back.col),
                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    legend.key = element_rect(fill = back.col)) +
              theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()) +
              theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank()) + facet_wrap(as.factor(data$MYconds))

          }
        }
      }
    }
################# Numeric
      if (is.numeric(data$clusters)) {
        if(length(MYConds) == 1) {
          myPLOT <- ggplot(data, aes(V6, y = V5,
                                     text = row.names(data), color = clusters)) +
            geom_point(size = cell.size, alpha = cell.transparency) +
            scale_colour_gradient(low = cell.colors[1], high = cell.colors[2], name="") +
            ggtitle(MyTitle) +
            theme(panel.background = element_rect(fill = back.col, colour = back.col),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  legend.key = element_rect(fill = back.col)) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())
        }
        if(length(MYConds) > 1) {
          myPLOT <- ggplot(data, aes(V6, y = V5,
                                     text = row.names(data), color = clusters)) +
            geom_point(size = cell.size, alpha = cell.transparency) +
            scale_colour_gradient(low = cell.colors[1], high = cell.colors[2], name="") +
            ggtitle(MyTitle) +
            theme(panel.background = element_rect(fill = back.col, colour = back.col),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  legend.key = element_rect(fill = back.col)) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()) + facet_wrap(as.factor(data$MYconds))

        }
      }
###################################
  # return
  if (interactive == TRUE) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}

