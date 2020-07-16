#' Plot nGenes, UMIs and percent mito
#'
#' This function takes an object of class iCellR and creates QC plot.
#' @param x An object of class iCellR.
#' @param plot.type Choose from "box.umi", "box.mito", "box.gene", "box.s.phase", "box.g2m.phase","all.in.one","three.in.one", "point.mito.umi", "point.gene.umi".
#' @param cell.color Choose a color for points in the plot.
#' @param cell.size A number for the size of the points in the plot, default = 1.
#' @param back.col Background color, default = "white"
#' @param box.color A color for the boxes in the "boxplot", default = "red".
#' @param box.line.col A color for the lines around the "boxplot", default = "green".
#' @param cell.transparency Color transparency for points in "scatterplot" and "boxplot", default = 0.5.
#' @param interactive If set to TRUE an interactive HTML file will be created, default = TRUE.
#' @param out.name If "interactive" is set to TRUE, the out put name for HTML, default = "plot".
#' @return An object of class iCellR.
#' @examples
#' stats.plot(demo.obj,
#'           plot.type = "three.in.one",
#'           out.name = "UMI-plot",
#'           interactive = FALSE,
#'           cell.color = "slategray3",
#'           cell.size = 1,
#'           cell.transparency = 0.5,
#'           box.color = "red",
#'           box.line.col = "green")
#' @import gridExtra
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly ggplotly layout plot_ly
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @importFrom methods new
#' @importFrom stats aggregate as.dendrogram cor cor.test dist hclust p.adjust prcomp quantile sd t.test
#' @importFrom utils capture.output packageVersion read.table write.table
#' @importFrom graphics legend par plot
#' @importFrom ggplot2 ggplot geom_segment geom_violin guide_colorbar guide_legend guides scale_color_discrete scale_colour_gradient scale_fill_gradient2 scale_x_continuous scale_y_continuous scale_y_discrete stat_summary coord_polar element_rect element_text element_blank facet_wrap scale_color_manual geom_hline geom_jitter geom_vline ylab xlab ggtitle theme_bw aes theme geom_bar geom_point geom_boxplot geom_errorbar position_dodge geom_tile geom_density geom_line
#' @export
stats.plot <- function (x = NULL,
                        plot.type = "three.in.one",
                        cell.color = "slategray3",
                        cell.size = 1,
                        cell.transparency = 0.5,
                        box.color = "red",
                        box.line.col = "green",
                        back.col = "white",
                        interactive = TRUE,
                        out.name = "plot")
{
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  DATA <- x@stats
  w = log1p(1:100)
  # get conditions
  do <- data.frame(do.call('rbind', strsplit(as.character(head(DATA$CellIds,1)),'_',fixed=TRUE)))
  do <- dim(do)[2]
if (do == 2) {
  col.legend <- data.frame(do.call('rbind', strsplit(as.character(DATA$CellIds),'_',fixed=TRUE)))[1]
  col.legend <- factor(as.matrix(col.legend))
} else {
  col.legend = "."
}
  # Box plots
  # mito
  mito.percent.plot <- ggplot(DATA,aes(y=mito.percent,x=col.legend)) +
    geom_jitter(height = 0,color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot( fill = box.color, col = "green", notch = FALSE, outlier.shape = NA, alpha = cell.transparency) +
    xlab("mito.percent") + ylab("percent of mito genes per cell") +
    stat_summary(fun=mean, geom="point", size=2, color="black") +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
    # nGenes
  nGenes.plot <- ggplot(DATA,aes(y=nGenes,x=col.legend)) +
    geom_jitter(height = 0,color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot( fill = box.color, col = box.line.col, notch = FALSE, outlier.shape = NA, alpha = cell.transparency) +
    xlab("nGenes") + ylab("number of genes per cell") +
    stat_summary(fun=mean, geom="point", size=2, color="black") +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
    # UMIs
  UMIsplot <- ggplot(DATA,aes(y=UMIs,x=col.legend)) +
    geom_jitter(height = 0,color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot( fill = box.color, col = box.line.col, notch = FALSE, outlier.shape = NA, alpha = cell.transparency) +
    xlab("UMIs") + ylab("number of UMIs per cell") +
    stat_summary(fun=mean, geom="point", size=2, color="black") +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
  # s.phase
  s.plot <- ggplot(DATA,aes(y=S.phase.probability,x=col.legend)) +
    geom_jitter(height = 0,color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot( fill = box.color, col = box.line.col, notch = FALSE, outlier.shape = NA, alpha = cell.transparency) +
    xlab("S phase") + ylab("S phase probability") +
    stat_summary(fun=mean, geom="point", size=2, color="black") +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
  # g2m.phase.probability
  g2m.plot <- ggplot(DATA,aes(y=g2m.phase.probability,x=col.legend)) +
    geom_jitter(height = 0,color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot( fill = box.color, col = box.line.col, notch = FALSE, outlier.shape = NA, alpha = cell.transparency) +
    xlab("G2 and M phase") + ylab("G2 and M phase probability") +
    stat_summary(fun=mean, geom="point", size=2, color="black") +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
# scatter plots
  if (col.legend[1] == ".") {
    Mito.UMIs <- ggplot(DATA,aes(y=mito.percent,x=UMIs,
                                 text = paste("UMIs =",UMIs,",",CellIds,sep=" "))) +
      geom_point(color = cell.color, size = cell.size, alpha = cell.transparency) +
      scale_x_continuous(trans = "log1p") +
      scale_color_discrete(name="") +
      theme(panel.background = element_rect(fill = back.col, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = back.col)) + theme_bw()
    #
    Genes.UMIs <- ggplot(DATA,aes(y=nGenes,x=UMIs,
                                  text = paste("nGenes =",nGenes,",",CellIds,sep=" "))) +
      geom_point(color = cell.color, size = cell.size, alpha = cell.transparency) +
      scale_x_continuous(trans = "log1p") +
      scale_y_continuous(trans = "log1p") +
      scale_color_discrete(name="") +
      theme(panel.background = element_rect(fill = back.col, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = back.col)) + theme_bw()
  } else {
    Mito.UMIs <- ggplot(DATA,aes(y=mito.percent,x=UMIs, col = col.legend,
                                 text = paste("UMIs =",UMIs,",",CellIds,sep=" "))) +
      geom_point( size = cell.size, alpha = cell.transparency) +
      scale_x_continuous(trans = "log1p") +
      scale_color_discrete(name="") +
      theme(panel.background = element_rect(fill = back.col, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = back.col)) + theme_bw()
    #
    Genes.UMIs <- ggplot(DATA,aes(y=nGenes,x=UMIs, col = col.legend,
                                  text = paste("nGenes =",nGenes,",",CellIds,sep=" "))) +
      geom_point(size = cell.size, alpha = cell.transparency) +
      scale_x_continuous(trans = "log1p") +
      scale_y_continuous(trans = "log1p") +
      scale_color_discrete(name="") +
      theme(panel.background = element_rect(fill = back.col, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = back.col)) + theme_bw()
  }
  # out puts
  if (plot.type == "point.mito.umi") {
    if (interactive == TRUE) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(Mito.UMIs),OUT.PUT)
    }
    else
    return(Mito.UMIs)
  }
  #
  if (plot.type == "point.gene.umi") {
    if (interactive == TRUE) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(Genes.UMIs),OUT.PUT)
    }
    else
    return(Genes.UMIs)
  }
  #
  if (plot.type == "box.umi") {
    if (interactive == TRUE) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(UMIsplot),OUT.PUT)
    }
    else
    return(UMIsplot)
  }
  #
  if (plot.type == "box.mito") {
    if (interactive == TRUE) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(mito.percent.plot),OUT.PUT)
    }
    else
    return(mito.percent.plot)
  }
  #
  if (plot.type == "box.gene") {
    if (interactive == TRUE) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(nGenes.plot),OUT.PUT)
    }
    else
    return(nGenes.plot)
  }
  #
  if (plot.type == "box.s.phase") {
    if (interactive == TRUE) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(s.plot),OUT.PUT)
    }
    else
      return(s.plot)
  }
  #
  if (plot.type == "box.g2m.phase") {
    if (interactive == TRUE) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(g2m.plot),OUT.PUT)
    }
    else
      return(g2m.plot)
  }
  #
  if (plot.type == "all.in.one") {
    if (interactive == TRUE) {
      message("for interactive mode use single plots (i.e. box.mito, box.gene, etc.) in plot.type")
    }
    return(grid.arrange(nGenes.plot,
                        UMIsplot,
                        mito.percent.plot,
                        s.plot,
                        g2m.plot,
                        ncol = 5))
  }
  #
  if (plot.type == "three.in.one") {
    if (interactive == TRUE) {
      message("for interactive mode use single plots (i.e. box.mito, box.gene, etc.) in plot.type")
    }
    return(grid.arrange(nGenes.plot,
                        UMIsplot,
                        mito.percent.plot,
                        ncol = 3))
  }
}
