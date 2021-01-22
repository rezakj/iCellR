#' Run diffusion map on PCA data (PHATE - Potential of Heat-Diffusion for Affinity-Based Transition Embedding)
#'
#' This function takes an object of class iCellR and runs diffusion map on PCA data.
#' @param x An object of class iCellR.
#' @param method diffusion map method, default = "phate".
#' @param dims PC dimentions to be used for UMAP analysis.
#' @param ndim int, optional, default: 2 number of dimensions in which the data will be embedded
#' @param k int, optional, default: 5 number of nearest neighbors on which to build kernel
#' @param alpha int, optional, default: 40 sets decay rate of kernel tails. If NULL, alpha decaying kernel is not used
#' @param n.landmark int, optional, default: 2000 number of landmarks to use in fast PHATE
#' @param gamma float, optional, default: 1 Informational distance constant between -1 and 1. gamma=1 gives the PHATE log potential, gamma=0 gives a square root potential.
#' @param t int, optional, default: 'auto' power to which the diffusion operator is powered sets the level of diffusion
#' @param knn.dist.method string, optional, default: 'euclidean'. recommended values: 'euclidean', 'cosine', 'precomputed' Any metric from scipy.spatial.distance can be used distance metric for building kNN graph. If 'precomputed', data should be an n_samples x n_samples distance or affinity matrix. Distance matrices are assumed to have zeros down the diagonal, while affinity matrices are assumed to have non-zero values down the diagonal. This is detected automatically using data[0,0]. You can override this detection with knn.dist.method='precomputed_distance' or knn.dist.method='precomputed_affinity'.
#' @param init phate object, optional object to use for initialization. Avoids recomputing intermediate steps if parameters are the same.
#' @param mds.method string, optional, default: 'metric' choose from 'classic', 'metric', and 'nonmetric' which MDS algorithm is used for dimensionality reduction
#' @param mds.dist.method string, optional, default: 'euclidean' recommended values: 'euclidean' and 'cosine'
#' @param t.max int, optional, default: 100. Maximum value of t to test for automatic t selection.
#' @param npca int, optional, default: 100 Number of principal components to use for calculating neighborhoods. For extremely large datasets, using n_pca < 20 allows neighborhoods to be calculated in log(n_samples) time.
#' @param plot.optimal.t boolean, optional, if TRUE, produce a plot showing the Von Neumann Entropy curve for automatic t selection.
#' @param verbose int or boolean, optional (default : 1) If TRUE or > 0, message verbose updates.
#' @param n.jobs int, optional (default: 1) The number of jobs to use for the computation. If -1 all CPUs are used. If 1 is given, no parallel computing code is used at all, which is useful for debugging. For n_jobs below -1, (n.cpus + 1 + n.jobs) are used. Thus for n_jobs = -2, all CPUs but one are used
#' @param seed int or NULL, random state (default: NULL)
#' @param potential.method Deprecated. For log potential, use gamma=1. For sqrt potential, use gamma=0.
#' @param use.alpha Deprecated To disable alpha decay, use alpha=NULL
#' @param n.svd Deprecated.
#' @param pca.method Deprecated.
#' @param g.kernel Deprecated.
#' @param diff.op Deprecated.
#' @param landmark.transitions Deprecated.
#' @param diff.op.t Deprecated.
#' @param dist.method Deprecated.
#' @return An object of class iCellR.
#' @export
run.diffusion.map <- function (x = NULL,
                      dims = 1:10,
                      method = "destiny",
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
  ###########
  ##########
  # https://github.com/lmcinnes/umap
  # get PCA data
  DATA <- x@pca.data
  DATA <- DATA[dims]
  #  2 dimention
  #  if (clust.dim == 2) {
  # TransPosed <- t(TopNormLogScale)
  if (method == "phate") {
    if(!"phateR" %in% (.packages())){
      stop("Please load phateR package: library(phateR)")
    }
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
  }
  if (method == "destiny") {
    if(!"destiny" %in% (.packages())){
      stop("Please load destiny package: library(destiny)")
    }
    DD <- DiffusionMap(DATA)
    HH <- as.data.frame(cbind(row.names(DATA),
                              as.numeric(DD$DC2),
                              as.numeric(DD$DC3),
                              as.numeric(DD$DC1)))
    colnames(HH) <- c("Ids","V1","V2","V3")
    row.names(HH) <- HH$Ids
    HH <- (HH)[,-1]
    DATA <- HH
    DATA$V1 <- as.numeric(scale(as.numeric(DATA$V1)))
    DATA$V2 <- as.numeric(scale(as.numeric(DATA$V2)))
    DATA$V3 <- as.numeric(scale(as.numeric(DATA$V3)))
  }
  attributes(x)$diffusion.data <- DATA
  # return
  return(x)
}
