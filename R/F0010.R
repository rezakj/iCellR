#' Make statistical information for each gene across all the cells (SD, mean, expression, etc.)
#'
#' This function takes an object of class iCellR and provides some statistical information for the genes.
#' @param x An object of class iCellR.
#' @param which.data Choose from "raw.data" or "main.data", default = "raw.data".
#' @param each.cond If TRUE each condition will be calculated, default = FALSE.
#' @return An object of class iCellR.
#' @examples
#' demo.obj <- gene.stats(demo.obj, which.data = "main.data")
#' head(demo.obj@gene.data)
#'
#' @export
gene.stats <- function (x = NULL,
                        which.data = "raw.data",
                        each.cond = FALSE) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
# get data
  if (which.data == "raw.data") {
    DATA <- x@raw.data
  }
  if (which.data == "main.data") {
    DATA <- x@main.data
  }
  # get conditions
  if (each.cond == TRUE) {
    do <- data.frame(do.call('rbind', strsplit(as.character(colnames(DATA)),'_',fixed=TRUE)))[1]
    Myconds <- as.character(as.matrix(unique(do)))
    if (length(Myconds) > 1) {
      for (i in Myconds) {
        ha <- subset(do, do == i)
        dim2 <- max(as.numeric(rownames(ha)))
        dim1 <- min(as.numeric(rownames(ha)))
        myDATA <- DATA[dim1:dim2]
        mymat = as.matrix(myDATA)
        SDs <- apply(mymat, 1, function(mymat) {sd(mymat)})
        Table <- list(row.names(myDATA),
                      as.numeric(rowSums(myDATA > 0)),
                      rep(dim(myDATA)[2], dim(myDATA)[1]),
                      (as.numeric(rowSums(myDATA > 0)) / dim(myDATA)[2])*100,
                      as.numeric(rowMeans(myDATA)),
                      as.numeric(SDs))
        names(Table) <- c("genes","numberOfCells","totalNumberOfCells","percentOfCells","meanExp","SDs")
        Table <- as.data.frame(Table)
        Table <- cbind(Table, condition = i)
        NameCol=paste("MyGeneStat",i,sep="_")
        eval(call("<-", as.name(NameCol), Table))
      }
    }
  }
  ###### do it for all too
  # calculate
  mymat = as.matrix(DATA)
  SDs <- apply(mymat, 1, function(mymat) {sd(mymat)})
  Table <- list(row.names(DATA),
                as.numeric(rowSums(DATA > 0)),
                rep(dim(DATA)[2], dim(DATA)[1]),
                (as.numeric(rowSums(DATA > 0)) / dim(DATA)[2])*100,
                as.numeric(rowMeans(DATA)),
                as.numeric(SDs))
  names(Table) <- c("genes","numberOfCells","totalNumberOfCells","percentOfCells","meanExp","SDs")
  Table <- as.data.frame(Table)
  MyGeneStat_all <- cbind(Table, condition = "all")
  ### merge
  filenames <- ls(pattern="MyGeneStat_")
  datalist <- mget(filenames)
  Table <- do.call(rbind.data.frame, datalist)
  row.names(Table) <- NULL
  attributes(x)$gene.data <- Table
  return(x)
}
