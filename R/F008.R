#' Normalize data
#'
#' This function takes an object of class iCellR and normalized the data based on "global.glsf", "ranked.glsf" or "spike.in" methods.
#' @param x An object of class iCellR.
#' @param norm.method Choose a normalization method, there are three option currently.
#' Choose from "global.glsf", "ranked.glsf","spike.in" or no.norm, default = "ranked.glsf".
#' @param top.rank If the method is set to "ranked.glsf", you need to set top number of genes sorted based on global base mean, default = 500.
#' @param spike.in.factors A numeric vector of spike-in values with the same cell id order as the main data.
#' @param rpm.factor If the norm.method is set to "rpm" the library sizes would be divided by this number, default = 1000 (higher numbers recomanded for bulk RNA-Seq).
#' @param rounding.digits integer indicating the number of decimal places (round) or significant digits (signif) to be used.
#' @param round.num Rounding of Numbers, default = FALSE.
#' @param ATAC.data If TURE, it would normalize ATAC-Seq data and not RNA-Seq, default = FALSE.
#' @param ATAC.filter If TURE, all the cells filtered in RNA-Seq will be filtered in ATAC-Seq. This needs to be done for both data to match,  default = TRUE.
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
                       rpm.factor = 1000,
                       rounding.digits = 3,
                       round.num = TRUE,
                       ATAC.data = FALSE,
                       ATAC.filter = TRUE) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
#######
# get data RNA
  DATA <- x@main.data
# get data ATAC
  if (ATAC.data == TRUE) {
    DATA <- x@atac.main
  }
#######
  if (norm.method == "global.glsf") {
    libSiz <- colSums(DATA)
    norm.facts <- as.numeric(libSiz) / mean(as.numeric(libSiz))
    dataMat <- as.matrix(DATA)
    normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  }
  if (norm.method == "ranked.glsf") {
    dataMat <- as.matrix(DATA)
    raw.data.order <- dataMat[ order(rowMeans(dataMat), decreasing = TRUE), ]
    topGenes = head(raw.data.order,top.rank)
    libSiz <- colSums(topGenes)
    norm.facts <- as.numeric(libSiz) / mean(as.numeric(libSiz))
    normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
  }
#######################
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
#####
  if (round.num == TRUE) {
    normalized <- round(normalized, digits = rounding.digits)
  }
######
  # get data ATAC
  if (ATAC.data == TRUE) {
    if (ATAC.filter == TRUE) {
      dat = normalized
      MyIDsToKeep <- colnames(x@main.data)
      normalized <- dat[ , which(names(dat) %in% MyIDsToKeep)]
      #### order colnames by the original data
      normalized <- as.matrix(normalized)
      normalized <- normalized[ , order(match(colnames(normalized),MyIDsToKeep))]
      normalized <- as.data.frame(normalized)
    }
    attributes(x)$atac.main <- normalized
  }
  # get data ATAC
  if (ATAC.data == FALSE) {
    attributes(x)$main.data <- normalized
    attributes(x)$norm.factors <- norm.facts
  }
############
  return(x)
}

