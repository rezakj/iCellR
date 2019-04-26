#' Load 10X data as data.frame
#'
#' This function takes 10X data files barcodes.tsv, genes.tsv and matrix.mtx and converts them to proper matrix file for iCellR.
#' @param dir.10x A directory that includes the 10X barcodes.tsv, genes.tsv and matrix.mtx files.
#' @param gene.name Should be either geneSymbol or ensembleID.
#' @return The data frame object
#' @examples
#' \dontrun{
#' load10x("/hg19", gene.name = "geneSymbol")
#' }
#' @import Matrix
#' @export
load10x <- function (dir.10x = NULL, gene.name = "geneSymbol") {
  if (!dir.exists(dir.10x)) {
    stop("Directory is not provided. Please provide a standard 10x matrix directory
         that includes; barcodes.tsv, genes.tsv and matrix.mtx files")
  }
  if (dir.exists(dir.10x)) {
  Standard10xInput <- paste(c("barcodes.tsv","genes.tsv","matrix.mtx"),collapse="")
  InputFiles <- paste(list.files(dir.10x),collapse="")
  }
  if (Standard10xInput != InputFiles) {
    stop("Provided directory does not have a standard 10x matrix directory
         that includes; barcodes.tsv, genes.tsv and matrix.mtx files")
  }
  if (Standard10xInput == InputFiles) {
    MTX10x <- list.files(dir.10x,full.names=T,"matrix.mtx")
    cell.barcodes <- list.files(dir.10x,full.names=T,"barcodes.tsv")
    gene.names.ids <- list.files(dir.10x,full.names=T,"genes.tsv")
    MTX10x <- readMM(MTX10x)
    cell.barcodes <- readLines(cell.barcodes)
  }
  if (gene.name == "geneSymbol") {
    gene.names.ids <- as.character(as.matrix(read.table(gene.names.ids,header=F)[2]))
  }
  if (gene.name == "ensembleID") {
    gene.names.ids <- as.character(as.matrix(read.table(gene.names.ids,header=F)[1]))
  }
  colnames(x = MTX10x) <- cell.barcodes
  rownames(x = MTX10x) <- gene.names.ids
  data.10x <- as.data.frame(as.matrix(MTX10x))
  row.names(data.10x) <- make.names(row.names(data.10x), unique=TRUE)
  row.names(data.10x) <- gsub("-",".",row.names(data.10x))
  return(data.10x)
}

