#' Gene-gene correlation.
#' This function helps to visulaize and calculate gene-gene correlations.
#' @param x An object of class iCellR.
#' @param data.type Choose from imputed and main, default = "imputed".
#' @param gene1 First gene name.
#' @param gene2 Second gene name.
#' @param conds Filter only one condition (only one), default is all conditions.
#' @param clusts Choose clusters to plot.
#' @param cell.size A numeric value for the size of the cells, default = 1.
#' @param cell.transparency A numeric value between 0 to 1, default = 0.5.
#' @param interactive If TRUE an html interactive file will be made, default = TRUE.
#' @param out.name Output name for html file if interactive = TRUE, default = "plot".
#' @return An object of class iCellR
#'
#' @export
gg.cor <- function (x = NULL,
                    data.type = "imputed",
                    gene1 = NULL,
                    gene2 = NULL,
                    conds = NULL,
                    clusts = NULL,
                    cell.size = 1,
                    cell.transparency = 0.5,
                    interactive = TRUE,
                    out.name = "plot") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ## get main data
  if (data.type == "main") {
    DATA <- x@main.data
  }
  if (data.type == "imputed") {
    DATA <- x@imputed.data
  }
  ###
  AllGenes = row.names(DATA)
  check1 = gene1 %in% AllGenes
  check2 = gene2 %in% AllGenes
  if (check1 == FALSE) {
    stop("gene 1 is not in the data")
  }
  if (check2 == FALSE) {
    stop("gene 2 is not in the data")
  }
  ####### Subset 1
  if (length(gene1) != 1) {
    stop("Please choose 1 gene")
  }
  if (length(gene2) != 1) {
    stop("Please choose 1 gene")
  }
  genes <- c(gene1, gene2)
  if (length(unique(genes)) != 2) {
    stop("Please choose 2 different genes for gene1 and gene2")
  }
  DATA <- subset(DATA, row.names(DATA) %in% genes)
  DATA <- as.data.frame(t(DATA))
#  head(DATA)
##### clusters
  col.legend <- as.data.frame(x@best.clust)
  DATA <- cbind(DATA, col.legend)
###### conditions
  ####### Subset 2
  MY.conds <- as.data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
  MY.conds <- as.character(as.matrix(MY.conds))
  DATA <- as.data.frame(cbind(DATA,MY.conds))
#############
  if (!is.null(conds)) {
    DATA <- subset(DATA, DATA$MY.conds %in% conds)
  }
  ####
  if (!is.null(clusts)) {
    DATA <- subset(DATA, DATA$clusters %in% clusts)
  }
########
  # cor
  MyCorelation <- as.numeric(cor(DATA[1],DATA[2]))
  MyCorelation <- round(MyCorelation, digits = 3)
  MyCorelation <-paste("correlation value:",MyCorelation)
  # pvalue
  PVal <- cor.test(as.numeric(as.matrix(DATA[1])),as.numeric(as.matrix(DATA[2])))$p.value
  if(PVal == 0) {
    PVal = "2.2e-16"
  }
  PVal <- paste("P value:", PVal)
  Sub <- paste(PVal, MyCorelation, sep= " | ")
  # plot
  Clusters <- factor(DATA$clusters)
  myPLOT <- ggplot(DATA, aes(x=log2(DATA[,2] +1 ), y=log2(DATA[,1] + 1), text = row.names(DATA), color = Clusters)) +
    geom_point(size = cell.size, alpha = cell.transparency) +
    xlab(paste("log2 expression (", gene1,")", sep="")) +
    ylab(paste("log2 expression (", gene2,")", sep="")) +
#    geom_abline(color= ab.col, linetype="dashed", size = ab.size) +
    ggtitle("gene gene correlation and cell gating", subtitle = Sub) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
#    scale_y_continuous(trans = "log1p") +
#    scale_x_continuous(trans = "log1p") +
    theme_bw()
##
  if (interactive == TRUE) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}
