#' Normalize ADT data.
#' This function takes data frame and Normalizes ADT data.
#' @param x An object of class iCellR.
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' gg.cor(my.obj, interactive = F, gene1 = "NKG7",gene2 = "GNLY", conds=c("WT"))
#' }
#'
#' @export
gg.cor <- function (x = NULL,
                    data.type = "imputed",
                    gene1 = NULL,
                    gene2 = NULL,
                    conds = NULL,
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
  check2 =gene2 %in% AllGenes
  if (check1 == F) {
    stop("gene 1 is not in the data")
  }
  if (check2 == F) {
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
  if (is.null(conds)) {
    DATA <- DATA
  } else {
    TAKE = grep(conds,row.names(DATA), value = T)
    DATA <- subset(DATA, row.names(DATA) %in% TAKE)
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
  myPLOT <- ggplot(DATA, aes(x=log2(DATA[,1] +1 ), y=log2(DATA[,2] + 1), text = row.names(DATA), color = Clusters)) +
    geom_point(size = cell.size, alpha = cell.transparency) +
    xlab(paste("log2 expression (", gene1,")", sep="")) +
    ylab(paste("log2 expression (", gene2,")", sep="")) +
#    geom_abline(color= ab.col, linetype="dashed", size = ab.size) +
    ggtitle("gene gene correlation and cell gating", subtitle = Sub) +
#    scale_y_continuous(trans = "log1p") +
#    scale_x_continuous(trans = "log1p") +
    theme_bw()
##
  if (interactive == T) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}
