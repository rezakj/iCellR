#' Calculate the number of UMIs and genes per cell and percentage of mitochondrial genes per cell and cell cycle genes.
#'
#' This function takes data frame and calculates the number of UMIs, genes per cell and percentage of mitochondrial genes per cell and cell cycle genes.
#' @param x A data frame containing gene counts for cells.
#' @param which.data Choose from raw data or main data, default = "raw.data".
#' @param mito.genes A character vector of mitochondrial  genes names , default is the genes starting with mt.
#' @param s.phase.genes A character vector of gene names for S phase, default = s.phase.
#' @param g2m.phase.genes A character vector of gene names for G2 and M phase, default = g2m.phase.
#' @return The data frame object
#' @examples
#' New.demo.obj <- qc.stats(demo.obj)
#' head(New.demo.obj@stats)
#' @export
qc.stats <- function (x = NULL,
                      which.data = "raw.data",
                      mito.genes = NULL,
                      s.phase.genes = s.phase,
                      g2m.phase.genes = g2m.phase) {
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
  # get UMIs
  UMIs <- colSums(DATA)
  # get nGENEs
  nGenes <- sapply(DATA, function(DATA) length(as.numeric(subset(DATA, DATA != 0))))
  # get mito gene names
  if (is.null(mito.genes)) {
    mito.genes <- grep(pattern = "^mt\\-", x = rownames(DATA), value = TRUE, ignore.case = TRUE)
    if ( length(mito.genes) == 0 ) {
      mito.genes <- grep(pattern = "^mt\\.", x = rownames(DATA), value = TRUE, ignore.case = TRUE)
    }
  }
  if (!is.null(mito.genes)) {
    mito.genes = mito.genes
  }
  # get mito percent
  mito <- subset(DATA,rownames(DATA) %in% mito.genes)
  if(dim(mito)[1] == 0){
    mitoSiz <- rep(0,dim(mito)[2])
    mito.percent <- mitoSiz
  }
  mitoSiz <- colSums(mito)
  mito.percent <- (mitoSiz/UMIs)
  # get cell cycle genes
  data <- row.names(DATA)
  # s
  s.phase.genes <- paste("^",s.phase.genes,"$", sep="")
  s.phase.genes <- paste(s.phase.genes,collapse="|")
  s.phase.genes <- grep(s.phase.genes, x = data, value = T, ignore.case = TRUE)
  s.phase.genes1 <- subset(DATA,rownames(DATA) %in% s.phase.genes)
  if(dim(s.phase.genes1)[1] != 0){
    s.phase.genes <- colSums(s.phase.genes1)
    s.phase.genes <- (s.phase.genes/UMIs)
  }
  if(dim(s.phase.genes1)[1] == 0){
    s.phase.genes <- rep(0,dim(s.phase.genes1)[2])
  }
  # g2m
  g2m.phase.genes <- paste("^",g2m.phase.genes,"$", sep="")
  g2m.phase.genes <- paste(g2m.phase.genes,collapse="|")
  g2m.phase.genes <- grep(g2m.phase.genes, x = data, value = T, ignore.case = TRUE)
  g2m.phase.genes1 <- subset(DATA,rownames(DATA) %in% g2m.phase.genes)
  if(dim(s.phase.genes1)[1] != 0){
    g2m.phase.genes <- colSums(g2m.phase.genes1)
    g2m.phase.genes <- (s.phase.genes/UMIs)
  }
  if(dim(g2m.phase.genes1)[1] == 0){
    g2m.phase.genes <- rep(0,dim(g2m.phase.genes1)[2])
  }
  # get all the data into DF
  QC <- list(colnames(DATA),
             as.numeric(nGenes),
             as.numeric(UMIs),
             as.numeric(mito.percent),
             as.numeric(s.phase.genes),
             as.numeric(g2m.phase.genes))
  names(QC) <- c("CellIds",
                 "nGenes",
                 "UMIs",
                 "mito.percent",
                 "S.phase.probability",
                 "g2m.phase.probability")
  STATS <- as.data.frame(QC)
  attributes(x)$stats <- STATS
  return(x)
}
