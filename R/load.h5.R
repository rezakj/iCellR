#' Load h5 data as data.frame
#'
#' This function reads hdf5 files.
#' @param filename path to the input (h5) file
#' @param feature.names row names to be feature names or ID numbers.
#' @param uniq.rows make row names unique.
#' @return The data frame object
#' @import hdf5r
#' @import Matrix
#' @import methods
#' @export
load.h5 <- function (filename, feature.names = TRUE, uniq.rows = TRUE)
{
  if (!file.exists(filename)) {
    stop("Input file not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = "r")
  genomes <- names(x = infile)
  output <- list()
  if (!infile$attr_exists("PYTABLES_FORMAT_VERSION")) {
    if (feature.names) {
      feature_slot <- "features/name"
    }
    else {
      feature_slot <- "features/id"
    }
  }
  else {
    if (feature.names) {
      feature_slot <- "gene_names"
    }
    else {
      feature_slot <- "genes"
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, "/data")]]
    indices <- infile[[paste0(genome, "/indices")]]
    indptr <- infile[[paste0(genome, "/indptr")]]
    shp <- infile[[paste0(genome, "/shape")]]
    features <- infile[[paste0(genome, "/", feature_slot)]][]
    barcodes <- infile[[paste0(genome, "/barcodes")]]
    sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[],
                               x = as.numeric(x = counts[]), dims = shp[], giveCsparse = FALSE)
    if (uniq.rows) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
    if (infile$exists(name = paste0(genome, "/features"))) {
      types <- infile[[paste0(genome, "/features/feature_type")]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(X = types.unique, FUN = function(x) {
          return(sparse.mat[which(x = types == x), ])
        }, simplify = FALSE, feature.names = TRUE)
      }
    }
    output[[genome]] <- sparse.mat
  }
######
#####
  infile$close_all()
  if (length(x = output) == 1) {
#    genome <- as.data.frame(as.matrix(genome))
    return(as.data.frame(as.matrix(output[[genome]])))
  }
  else {
    return(as.data.frame(as.matrix(output)))
  }
}
