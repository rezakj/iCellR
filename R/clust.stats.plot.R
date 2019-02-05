#' Plot nGenes, UMIs and perecent mito
#'
#' This function takes an object of class iCellR and creats QC plot.
#' @param x An object of class iCellR.
#' @param plot.type Choose from "box.umi", "box.mito", "box.gene", defult = "box.mito".
#' @param cell.color Choose a color for points in the plot.
#' @param cell.size A number for the size of the points in the plot, defult = 1.
#' @param box.color A color for the boxes in the "boxplot", defult = "red".
#' @param box.line.col A color for the lines around the "boxplot", defult = "green".
#' @param notch Notch the box plots, defult = F.
#' @param cell.transparency Color transparency for points in "scatterplot" and "boxplot", defult = 0.5.
#' @param interactive If set to TRUE an intractive HTML file will be created, defult = TRUE.
#' @param out.name If "interactive" is set to TRUE, the out put name for HTML, defult = "plot".
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' clust.stats.plot(my.obj, plot.type = "box.mito", interactive = F, out.name = "box.mito.clusters")
#' }
#' @export
clust.stats.plot <- function (x = NULL,
                        plot.type = "box.mito",
                        cell.color = "slategray3",
                        cell.size = 1,
                        cell.transparency = 0.5,
                        box.color = "red",
                        box.line.col = "green",
                        back.col = "white",
                        notch = F,
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
  DATA <- merge(DATA, MyClusts, by = "row.names", all.x=F, all.y=T)
  # plot
  # mito
  mito.percent.plot <- ggplot(DATA,aes(y=mito.percent, x=as.factor(clusters))) +
    geom_jitter(color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot(fill = box.color, col = "green", notch = notch, outlier.shape = NA, alpha = cell.transparency) +
    xlab("clusters") + ylab("percent of mito genes per cell") +
    stat_summary(fun.y=mean, geom="point", size=2, color="black") +
    theme_bw()
  # nGenes
  nGenes.plot <- ggplot(DATA,aes(y=nGenes,x=as.factor(clusters))) +
    geom_jitter(color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot(fill = box.color, col = box.line.col, notch = notch, outlier.shape = NA, alpha = cell.transparency) +
    xlab("clusters") + ylab("number of genes per cell") +
    stat_summary(fun.y=mean, geom="point", size=2, color="black") +
    theme_bw()
  # UMIs
  UMIsplot <- ggplot(DATA,aes(y=UMIs,x=as.factor(clusters))) +
    geom_jitter(color = cell.color, size = cell.size, alpha = cell.transparency) +
    geom_violin(trim=FALSE, col = "black", alpha = cell.transparency) +
    geom_boxplot(fill = box.color, col = box.line.col, notch = notch, outlier.shape = NA, alpha = cell.transparency) +
    xlab("clusters") + ylab("number of UMIs per cell") +
    stat_summary(fun.y=mean, geom="point", size=2, color="black") +
    theme_bw()
# return
  if (plot.type == "box.umi") {
    if (interactive == T) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(Mito.UMIs),OUT.PUT)
    }
    else
      return(UMIsplot)
  }
  #
  if (plot.type == "box.mito") {
    if (interactive == T) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(Mito.UMIs),OUT.PUT)
    }
    else
      return(mito.percent.plot)
  }
  #
  if (plot.type == "box.gene") {
    if (interactive == T) {
      OUT.PUT <- paste(out.name, ".html", sep="")
      htmlwidgets::saveWidget(ggplotly(Mito.UMIs),OUT.PUT)
    }
    else
      return(nGenes.plot)
  }
}
