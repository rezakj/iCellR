#' Load 10X data as data.frame
#'
#' This function takes 10X data files barcodes.tsv, genes.tsv and matrix.mtx and converts them to proper matrix file for iCellR.
#' @param dir.10x A directory that includes the 10X barcodes.tsv, genes.tsv and matrix.mtx files.
#' @param gene.name Gene names or ids column number, default = 2.
#' @return The data frame object
#' @examples
#' my.data <- load10x(system.file("extdata", "filtered_gene_bc_matrices", package = "iCellR"))
#'
#' # See first few rows and columns
#' head(my.data)[1:5]
#'
#' @import Matrix
#' @import knitr
#' @export
load10x <- function (dir.10x = NULL, gene.name = 2) {
  if (!dir.exists(dir.10x)) {
    stop("Directory is not provided. Please provide a standard 10x matrix directory
         that includes; barcodes.tsv, genes.tsv/features.tsv and matrix.mtx files (data could be zipped or unzipped)")
  }
    MTX10x <- list.files(dir.10x,full.names=T,pattern =c("matrix"))
    cell.barcodes <- list.files(dir.10x,full.names=T,pattern =c("barcodes"))
    gene.names.ids1 <- list.files(dir.10x,full.names=T,pattern =c("genes"))
    gene.names.ids2 <- list.files(dir.10x,full.names=T,pattern =c("features"))
    gene.names.ids <- c(gene.names.ids1,gene.names.ids2)
    MTX10x <- readMM(MTX10x)
    cell.barcodes <- readLines(cell.barcodes)
    cell.barcodes <- gsub("-",".",cell.barcodes)
  ## get gene name or ID
    gene.names.ids <- as.character(as.matrix(read.table(gene.names.ids,header=F,sep="\t")[gene.name]))
    gene.names.ids <- gsub("-",".",gene.names.ids)
    ##
  colnames(x = MTX10x) <- cell.barcodes
  rownames(x = MTX10x) <- gene.names.ids
  data.10x <- as.data.frame(as.matrix(MTX10x))
  row.names(data.10x) <- make.names(row.names(data.10x), unique=TRUE)
  return(data.10x)
}

