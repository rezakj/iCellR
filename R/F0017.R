#' Run tSNE on PCA Data. Barnes-Hut implementation of t-Distributed Stochastic Neighbor Embedding
#'
#' This function takes an object of class iCellR and runs tSNE on PCA data. Wrapper for the C++ implementation of Barnes-Hut t-Distributed Stochastic Neighbor Embedding. t-SNE is a method for constructing a low dimensional embedding of high-dimensional data, distances or similarities. Exact t-SNE can be computed by setting theta=0.0.
#' @param x An object of class iCellR.
#' @param dims PC dimentions to be used for tSNE analysis.
#' @param add.3d Add 3D tSNE as well, default = TRUE.
#' @param my.seed seed number, default = 0.
#' @param initial_dims integer; the number of dimensions that should be retained in the initial PCA step (default: 50)
#' @param perplexity numeric; Perplexity parameter
#' @param theta numeric; Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
#' @param check_duplicates logical; Checks whether duplicates are present. It is best to make sure there are no duplicates present and set this option to FALSE, especially for large datasets (default: TRUE)
#' @param pca logical; Whether an initial PCA step should be performed (default: TRUE)
#' @param max_iter integer; Number of iterations (default: 1000)
#' @param verbose logical; Whether progress updates should be messageed (default: FALSE)
#' @param is_distance logical; Indicate whether X is a distance matrix (experimental, default: FALSE)
#' @param Y_init matrix; Initial locations of the objects. If NULL, random initialization will be used (default: NULL). Note that when using this, the initial stage with exaggerated perplexity values and a larger momentum term will be skipped.
#' @param pca_center logical; Should data be centered before pca is applied? (default: TRUE)
#' @param pca_scale logical; Should data be scaled before pca is applied? (default: FALSE)
#' @param stop_lying_iter integer; Iteration after which the perplexities are no longer exaggerated (default: 250, except when Y_init is used, then 0)
#' @param mom_switch_iter integer; Iteration after which the final momentum is used (default: 250, except when Y_init is used, then 0)
#' @param momentum numeric; Momentum used in the first part of the optimization (default: 0.5)
#' @param final_momentum numeric; Momentum used in the final part of the optimization (default: 0.8)
#' @param eta numeric; Learning rate (default: 200.0)
#' @param exaggeration_factor numeric; Exaggeration factor used to multiply the P matrix in the first part of the optimization (default: 12.0)
#' @return An object of class iCellR.
#' @examples
#' demo.obj <- run.pc.tsne(demo.obj, dims = 1:10,perplexity = 20)
#'
#' head(demo.obj@pca.data)[1:5]
#'
#' @import Rtsne
#' @export
run.pc.tsne <- function (x = NULL,
                      dims = 1:10, my.seed = 0,add.3d = TRUE,
                      initial_dims = 50, perplexity = 30,
                      theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,
                      verbose = FALSE, is_distance = FALSE, Y_init = NULL,
                      pca_center = TRUE, pca_scale = FALSE,
                      stop_lying_iter = ifelse(is.null(Y_init), 250L, 0L),
                      mom_switch_iter = ifelse(is.null(Y_init), 250L, 0L), momentum = 0.5,
                      final_momentum = 0.8, eta = 200, exaggeration_factor = 12) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
# get PCA data
  set.seed(my.seed)
  DATA <- x@pca.data
  TransPosed <- DATA[dims]
  #  2 dimention
  #  if (clust.dim == 2) {
  # TransPosed <- t(TopNormLogScale)
  tsne <- Rtsne(TransPosed, dims = 2,
                initial_dims = initial_dims, perplexity = perplexity,
                theta = theta, check_duplicates = check_duplicates, pca = pca, max_iter = max_iter,
                verbose = verbose, is_distance = is_distance, Y_init = Y_init,
                pca_center = pca_center, pca_scale = pca_scale,
                stop_lying_iter = stop_lying_iter,
                mom_switch_iter = mom_switch_iter, momentum = momentum,
                final_momentum = final_momentum, eta = eta, exaggeration_factor = exaggeration_factor)
  tsne.data = as.data.frame(tsne$Y)
  tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
  rownames(tsne.data) <- tsne.data$cells
  tsne.data <- tsne.data[,-1]
  attributes(x)$tsne.data <- tsne.data
  #  }
  # choose 3 demention
  # tSNE
  # TransPosed <- t(TopNormLogScale)
  if (add.3d == TRUE) {
  tsne <- Rtsne(TransPosed, dims = 3,
                initial_dims = initial_dims, perplexity = perplexity,
                theta = theta, check_duplicates = check_duplicates, pca = pca, max_iter = max_iter,
                verbose = verbose, is_distance = is_distance, Y_init = Y_init,
                pca_center = pca_center, pca_scale = pca_scale,
                stop_lying_iter = stop_lying_iter,
                mom_switch_iter = mom_switch_iter, momentum = momentum,
                final_momentum = final_momentum, eta = eta, exaggeration_factor = exaggeration_factor)
  tsne.data = as.data.frame(tsne$Y)
  tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
  rownames(tsne.data) <- tsne.data$cells
  tsne.data <- tsne.data[,-1]
  attributes(x)$tsne.data.3d <- tsne.data
  }
  return(x)
}
