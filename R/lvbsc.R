#' Low Variance Batch Spacing Correction (LVBSC)
#'
#' This function takes an object of class iCellR and finds the top genes with low variance across the samples to make a model for batch spacing correction.
#' @param x An object of class iCellR.
#' @param dispersion.limit.max A number for taking the genes that have dispersion below this number, defult = 1.5.
#' @param base.mean.rank A number taking the top genes ranked by base mean, defult = 500.
#' @param non.sig.col Color for the genes not used for the model, defult = "darkgray".
#' @param right.sig.col Color for the genes above the dispersion limit, defult = "chartreuse3".
#' @param left.sig.col Color for the genes above the rank limit, defult = "cadetblue3".
#' @param disp.line.col Color of the line for dispersion limit, defult = "black".
#' @param rank.line.col Color of the line for rank limit, defult = "red".
#' @param cell.size A number for the size of the points in the plot, defult = 1.75.
#' @param no.mito.model If set to TRUE, mitochondrial genes would be excluded from the gene list made for clustering, defult = TRUE.
#' @param mark.mito Mark mitochondrial genes in the plot, defult = TRUE.
#' @param cell.transparency Color transparency for the points in the plot, defult = 0.5.
#' @param interactive If set to TRUE an intractive HTML file will be created, defult = TRUE.
#' @param out.name If "interactive" is set to TRUE, the out put name for HTML, defult = "plot".
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' lvbsc(my.obj,
#'                dispersion.limit.max = 1.5,
#'                base.mean.rank = 500,
#'                no.mito.model = T,
#'                mark.mito = T,
#'                interactive = T,
#'                out.name = "gene.model")
#'
#' lvbsc(my.obj,
#'              dispersion.limit.max = 1.5,
#'              base.mean.rank = 500)
#' }
#'
#' @import ggrepel
#' @export
lvbsc <- function (x = NULL,
                             dispersion.limit.max = 1.5,
                             base.mean.rank = 500) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  data <- x@gene.data
  #  x = c(1:100)
  #  y = log1p(1:100)
  # variables
  data <- subset(data, SDs < dispersion.limit.max)
  data <- data[order(data$meanExp, decreasing = T),]
  # get mito genes
  mito.genes <- grep(pattern = "^mt\\-", x = data$genes, value = TRUE, ignore.case = TRUE)
  if ( length(mito.genes) == 0 ) {
    mito.genes <- grep(pattern = "^mt\\.", x = data$genes, value = TRUE, ignore.case = TRUE)
  }
  my.clust.genes = head(data,base.mean.rank)
  my.clust.genes <- subset(my.clust.genes, !(genes %in% mito.genes))
  my.clust.genes <- (my.clust.genes$genes)
  write.table((my.clust.genes),file="lvbsc.txt", row.names =F, quote = F, col.names = F)
}
