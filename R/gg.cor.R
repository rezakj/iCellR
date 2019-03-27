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
                    ab.size = 1,
                    ab.col = "darkred",
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
  ################################### conditions
    col.legend <- data.frame(do.call('rbind', strsplit(as.character(row.names(DATA)),'_',fixed=TRUE)))[1]
    colnames(col.legend) <- "conditions"
####
  DATA <- cbind(DATA, col.legend)
  ####### Subset 2
  if (is.null(conds)) {
    DATA <- DATA
  } else {
    DATA <- subset(DATA, DATA[,3] %in% conds)
  }
########
#  LowessLine <- as.data.frame(lowess(DATA[1:2]))
#  DATA <- cbind(DATA, LowessLine)
#  X <- DATA$x
#  Y <- DATA$y
  # cor
  MyCorelation <- as.numeric(cor(DATA[1],DATA[2]))
  MyCorelation <- round(MyCorelation, digits = 3)
  MyCorelation <-paste("correlation value:",MyCorelation)
  # pvalue
  PVal <- cor.test(as.numeric(as.matrix(DATA[1])),as.numeric(as.matrix(DATA[2])))$p.value
  PVal <- paste("P value:", PVal)
  Sub <- paste(PVal, MyCorelation, sep= " | ")
  # plot
  myPLOT <- ggplot(DATA, aes(x=DATA[,1], y=DATA[,2], text = row.names(DATA), color = conditions)) +
    geom_point(size = cell.size, alpha = cell.transparency) +
    xlab(gene1) +
    ylab(gene2) +
    geom_abline(color= ab.col, linetype="dashed", size = ab.size) +
    ggtitle("gene gene correlation", subtitle = Sub) +
    theme_bw()
##
  if (interactive == T) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}
