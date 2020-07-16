#' Plotting tSNE, PCA, UMAP, Diffmap and other dim reductions
#'
#' This function takes an object of class iCellR and creates QC plot.
#' @param x An object of class iCellR.
#' @param plot.type Choose from "bar.cc", "pie.cc" , box.umi", "box.mito", "box.gene", default = "box.mito".
#' @param cell.color Choose a color for points in the plot.
#' @param cell.size A number for the size of the points in the plot, default = 1.
#' @param box.color A color for the boxes in the "boxplot", default = "red".
#' @param box.line.col A color for the lines around the "boxplot", default = "green".
#' @param notch Notch the box plots, default = FALSE.
#' @param back.col Background color, default = "white"
#' @param cell.transparency Color transparency for points in "scatterplot" and "boxplot", default = 0.5.
#' @param interactive If set to TRUE an interactive HTML file will be created, default = TRUE.
#' @param out.name If "interactive" is set to TRUE, the out put name for HTML, default = "plot".
#' @param conds.to.plot Choose the conditions you want to see in the plot, default = NULL (all conditions).
#' @return An object of class iCellR.
#' @examples
#' clust.stats.plot(demo.obj,
#'                    plot.type = "box.mito",
#'                    interactive = FALSE,
#'                    out.name = "box.mito.clusters")
#' @export
clust.stats.plot <- function (x = NULL,
                        plot.type = "box.mito",
                        conds.to.plot = NULL,
                        cell.color = "slategray3",
                        cell.size = 1,
                        cell.transparency = 0.5,
                        box.color = "red",
                        box.line.col = "green",
                        back.col = "white",
                        notch = FALSE,
                        interactive = TRUE,
                        out.name = "plot")
{
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  # get stats data for all cells
  DATA <- x@stats
  row.names(DATA) <- DATA$CellIds
  # get cluster data for all cells
  MyClusts <- x@best.clust
  # merge
  DATA <- merge(DATA, MyClusts, by = "row.names", all.x=FALSE, all.y=TRUE)
##### conditions
  Cells <- as.character(DATA$CellIds)
  MYConds <- data.frame(do.call('rbind', strsplit(as.character(Cells),'_',fixed=TRUE)))[1]
  colnames(MYConds) <- "conditions"
  DATA <- cbind(DATA,MYConds)
  if (!is.null(conds.to.plot)) {
    DATA <- subset(DATA, DATA$conditions %in% conds.to.plot)
  }
  ######  # plot
# cell cycle
  # bar
  COUNTS <- c(1:length(DATA$Phase))
  myBP <- ggplot(DATA,aes(y=COUNTS,
                          x=as.factor(clusters), fill = Phase)) +
    geom_bar(stat = "identity") + theme_bw() +
    theme(axis.text.x=element_text(angle=90)) +
    ylab("Cell number ratio")
  # pie
  myPIE <- ggplot(DATA,aes(y=clusters, x="", fill = Phase)) +
    geom_bar(stat = "identity", position = "fill") + theme_bw() + facet_wrap(~ clusters) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + coord_polar(theta="y")
  # mito
  mito.percent.plot <- ggplot(DATA,aes(y=mito.percent, x=as.factor(clusters))) +
    geom_jitter(height = 0,color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot(fill = box.color, col = "green", notch = notch, outlier.shape = NA, alpha = cell.transparency) +
    xlab("clusters") + ylab("percent of mito genes per cell") +
    stat_summary(fun=mean, geom="point", size=2, color="black") +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
  # nGenes
  nGenes.plot <- ggplot(DATA,aes(y=nGenes,x=as.factor(clusters))) +
    geom_jitter(height = 0,color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot(fill = box.color, col = box.line.col, notch = notch, outlier.shape = NA, alpha = cell.transparency) +
    xlab("clusters") + ylab("number of genes per cell") +
    stat_summary(fun=mean, geom="point", size=2, color="black") +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
  # UMIs
  UMIsplot <- ggplot(DATA,aes(y=UMIs,x=as.factor(clusters))) +
    geom_jitter(height = 0,color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot(fill = box.color, col = box.line.col, notch = notch, outlier.shape = NA, alpha = cell.transparency) +
    xlab("clusters") + ylab("number of UMIs per cell") +
    stat_summary(fun=mean, geom="point", size=2, color="black") +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
# return
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
  if (plot.type == "pie.cc") {
    if (interactive == TRUE) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(myPIE),OUT.PUT)
    }
    else
      return(myPIE)
  }
  #
  if (plot.type == "bar.cc") {
    if (interactive == TRUE) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(myBP),OUT.PUT)
    }
    else
      return(myBP)
  }
}
