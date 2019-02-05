#' Run PHATE on PCA data (PHATE - Potential of Heat-Diffusion for Affinity-Based Transition Embedding)
#'
#' This function takes an object of class iCellR and runs PHATE on PCA data.
#' @param x An object of class iCellR.
#' @param dims PC dimentions to be used for UMAP analysis.
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- run.diffusion.map(my.obj, dims = 1:10, method = "phate")
#' }
#' @export
run.diffusion.map <- function (x = NULL,
                      dims = 1:10,
                      method = "phate",
                      ndim = 3,
                      k = 5, alpha = 40, n.landmark = 2000,
                      gamma = 1, t = "auto", knn.dist.method = "euclidean",
                      init = NULL, mds.method = "metric", mds.dist.method = "euclidean",
                      t.max = 100, npca = 100, plot.optimal.t = FALSE, verbose = 1,
                      n.jobs = 1, seed = NULL, potential.method = NULL,
                      use.alpha = NULL, n.svd = NULL, pca.method = NULL,
                      g.kernel = NULL, diff.op = NULL, landmark.transitions = NULL,
                      diff.op.t = NULL, dist.method = NULL) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  # https://github.com/lmcinnes/umap
  # get PCA data
  require(phateR)
  DATA <- x@pca.data
  DATA <- DATA[dims]
  #  2 dimention
  #  if (clust.dim == 2) {
  # TransPosed <- t(TopNormLogScale)
  if (method == "phate") {
  DD <- phate(DATA, ndim = ndim, k = k, alpha = alpha, n.landmark = n.landmark,
              gamma = gamma, t = t, knn.dist.method = knn.dist.method,
              init = init, mds.method = mds.method, mds.dist.method = mds.dist.method,
              t.max = t.max, npca = npca, plot.optimal.t = plot.optimal.t, verbose = verbose,
              n.jobs = n.jobs, seed = seed, potential.method = potential.method,
              use.alpha = use.alpha, n.svd = n.svd, pca.method = pca.method,
              g.kernel = g.kernel, diff.op = diff.op, landmark.transitions = landmark.transitions,
              diff.op.t = diff.op.t, dist.method = dist.method)
  DATA <- as.data.frame(DD$embedding)
#  DATA <- as.data.frame(scale(DD$embedding))
  attributes(x)$diffusion.data <- DATA
  # return
  return(x)
  }
}
