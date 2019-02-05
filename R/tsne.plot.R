#' Plot nGenes, UMIs and perecent mito
#'
#' This function takes an object of class iCellR and creats plot.
#' @param x An object of class iCellR.
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' tsne.plot(my.obj)
#' }
#' @export
tsne.plot <- function (x = NULL, cell.size = 1) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  DATA <- x@tsne.data
  clusters <- DATA$cl_hierarchical
  myPLOT <- ggplot(DATA, aes(x = V1, y = V2,
                             text = cells ,
                             color = clusters)) +
    geom_point(size = cell.size, alpha = 0.5) +
    xlab("Dim1") +
    ylab("Dim2") +
    ggtitle("t-SNE plot") + theme_bw()
  return(myPLOT)
}
