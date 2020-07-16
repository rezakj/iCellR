#' Normalize data
#'
#' This function takes an object of class iCellR and normalized the data based on "global.glsf", "ranked.glsf" or "spike.in" methods.
#' @param x An object of class iCellR.
#' @param norm.method Choose a normalization method, there are three option currently.
#' Choose from "global.glsf", "ranked.glsf","spike.in" or no.norm, default = "ranked.glsf".
#' @param top.rank If the method is set to "ranked.glsf", you need to set top number of genes sorted based on global base mean, default = 500.
#' @param spike.in.factors A numeric vector of spike-in values with the same cell id order as the main data.
#' @param rpm.factor If the norm.method is set to "rpm" the library sizes would be divided by this number, default = 1000 (higher numbers recomanded for bulk RNA-Seq).
#' @return An object of class iCellR.
#' @examples
#'
#' demo.obj <- norm.data(demo.obj, norm.method = "ranked.glsf", top.rank = 500)
#'
#' @export
norm.data <- function (x = NULL,
                       norm.method = "ranked.glsf",
                       top.rank = 500,
                       spike.in.factors = NULL,
                       rpm.factor = 1000) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  DATA <- x@main.data
  if (norm.method == "global.glsf") {
    libSiz <- colSums(DATA)
    norm.facts <- as.numeric(libSiz) / mean(as.numeric(libSiz))
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  }
  if (norm.method == "ranked.glsf") {
    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = T), ]
    topGenes = head(raw.data.order,top.rank)
    libSiz <- colSums(topGenes)
    norm.facts <- as.numeric(libSiz) / mean(as.numeric(libSiz))
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  }
#  if (norm.method == "deseq") {
#    require(DESeq)
#    CondAnum <- length(colnames(DATA)) - 1
#    conds <- factor( c( rep("CondA", CondAnum) , rep("CondB", 1)))
#    cds <- newCountDataSet(DATA, conds )
#    cds <- estimateSizeFactors( cds )
#    norm.facts <- as.numeric(sizeFactors(cds))
#    normalized <- as.data.frame(counts(cds,normalized=TRUE))
#  }
#  if (norm.method == "ranked.deseq") {
#    require(DESeq)
#    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = T), ]
#    CondAnum <- length(colnames(DATA)) - 1
#    conds <- factor( c( rep("CondA", CondAnum) , rep("CondB", 1)))
#    cds <- newCountDataSet(raw.data.order, conds )
#    cds <- estimateSizeFactors( cds )
#    norm.facts <- as.numeric(sizeFactors(cds))
#    cds1 <- newCountDataSet(DATA, conds)
#    sizeFactors(cds1) = sizeFactors(cds)
#    rm("cds")
#    normalized <- as.data.frame(counts(cds1, normalized=TRUE))
#    rm("cds1")
#  }
  if (norm.method == "spike.in") {
    norm.facts <- spike.in.factors
    norm.facts <- as.numeric(norm.facts)
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  }
  if (norm.method == "no.norm") {
    norm.facts <- colnames(DATA)
    norm.facts <- norm.facts == 1
    norm.facts[ norm.facts == "FALSE" ] <- 1
    normalized <- DATA
  }
  if (norm.method == "rpm") {
    libSiz <- colSums(DATA)
    norm.facts <- as.numeric(libSiz) / rpm.factor
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  }
#  SizeFactors <- as.numeric(norm.facts)
#  names(SizeFactors) <- c(colnames(DATA))
#  SizeFactors <- as.data.frame(SizeFactors)
#  normalized[is.na(normalized)] <- 0
  normalized <- round(normalized, digits = 3)
  attributes(x)$main.data <- normalized
  attributes(x)$norm.factors <- norm.facts
  return(x)
}

