#' Calculate cluster and conditions frequencies
#'
#' This function takes an object of class iCellR and calculates cluster and conditions frequencies.
#' @param x An object of class iCellR.
#' @param plot.type Choose from pie or bar, defult = pie.
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' clust.cond.info(my.obj, plot.type = "pie")
#' clust.cond.info(my.obj, plot.type = "bar")
#' }
#' @export
clust.cond.info <- function (x = NULL, plot.type = "pie" ) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ###################
  if (my.obj@data.conditions == "Data conditions: no conditions/single sample") {
    stop("This function only runs when there are conditions/more than one sample")
  }
  ###################
  DATA <- (x@best.clust)
  Conds <- (as.data.frame(do.call("rbind", strsplit(row.names(DATA), "_")))[1])
  clusts <- (as.data.frame(DATA$clusters))
  cond.clust <- cbind(Conds, clusts)
  colnames(cond.clust) <- c("conditions","clusters")
  DATA <- as.data.frame(table(cond.clust))
# as.data.frame(table(Conds))
  # bar
  myBP <- ggplot(DATA,aes(y=Freq, x=clusters, fill = conditions)) +
   geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x=element_text(angle=90))
  # pie
  myPIE <- ggplot(DATA,aes(y=Freq, x="", fill = conditions)) +
    geom_bar(stat = "identity", position = "fill") + theme_bw() + facet_wrap(~ clusters) +
      theme(axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) + coord_polar(theta="y")
  #############
  write.table((DATA), file="clust_cond_freq_info.txt", sep="\t", row.names =F)
  print("clust_cond_freq_info.txt file has beed generated.")
  if (plot.type == "bar") {
    return(myBP)
  }
  if (plot.type == "pie") {
    return(myPIE)
  }
#
}
