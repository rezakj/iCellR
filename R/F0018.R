#' Run UMAP on PCA Data (Computes a manifold approximation and projection)
#'
#' This function takes an object of class iCellR and runs UMAP on PCA data.
#' @param x An object of class iCellR.
#' @param dims PC dimentions to be used for UMAP analysis.
#' @param my.seed seed number, default = 0.
#' @param n_neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100.
#' @param n_components The dimension of the space to embed into. This defaults to 2 to provide easy visualization, but can reasonably be set to any integer value in the range 2 to 100.
#' @param metric Type of distance metric to use to find nearest neighbors. "euclidean" (the default)
#' @param n_epochs Number of epochs to use during the optimization of the embedded coordinates. By default, this value is set to 500 for datasets containing 10,000 vertices or less, and 200 otherwise. If n_epochs = 0, then coordinates determined by "init" will be returned.
#' @param learning_rate Initial learning rate used in optimization of the coordinates.
#' @param scale Scaling to apply to X if it is a data frame or matrix: "none" or FALSE or NULL No scaling. "Z" or "scale" or TRUE Scale each column to zero mean and variance 1.
#' @param init Type of initialization for the coordinates.
#' @param init_sdev If non-NULL, scales each dimension of the initialized coordinates (including any user-supplied matrix) to this standard deviation. By default no scaling is carried out, except when init = "spca", in which case the value is 0.0001. Scaling the input may help if the unscaled versions result in initial coordinates with large inter-point distances or outliers. This usually results in small gradients during optimization and very little progress being made to the layout. Shrinking the initial embedding by rescaling can help under these circumstances. Scaling the result of init = "pca" is usually recommended and init = "spca" as an alias for init = "pca", init_sdev = 1e-4 but for the spectral initializations the scaled versions usually aren't necessary unless you are using a large value of n_neighbors (e.g. n_neighbors = 150 or higher).
#' @param spread The effective scale of embedded points. In combination with min_dist, this determines how clustered/clumped the embedded points are.
#' @param min_dist The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. The value should be set relative to the spread value, which determines the scale at which embedded points will be spread out.
#' @param set_op_mix_ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets. Both fuzzy set operations use the product t-norm. The value of this parameter should be between 0.0 and 1.0; a value of 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy intersection.
#' @param local_connectivity The local connectivity required – i.e. the number of nearest neighbors that should be assumed to be connected at a local level. The higher this value the more connected the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.
#' @param bandwidth The effective bandwidth of the kernel if we view the algorithm as similar to Laplacian Eigenmaps. Larger values induce more connectivity and a more global view of the data, smaller values concentrate more locally.
#' @param repulsion_strength Weighting applied to negative samples in low dimensional embedding optimization. Values higher than one will result in greater weight being given to negative samples.
#' @param negative_sample_rate The number of negative edge/1-simplex samples to use per positive edge/1-simplex sample in optimizing the low dimensional embedding.
#' @param a More specific parameters controlling the embedding. If NULL these values are set automatically as determined by min_dist and spread.
#' @param b More specific parameters controlling the embedding. If NULL these values are set automatically as determined by min_dist and spread.
#' @param nn_method Method for finding nearest neighbors. Options are: "fnn". Use exact nearest neighbors via the FNN package. "annoy" Use approximate nearest neighbors via the RcppAnnoy package. "idx". A n_vertices x n_neighbors matrix containing the integer indexes of the nearest neighbors in X. Each vertex is considered to be its own nearest neighbor, i.e. idx[, 1] == 1:n_vertices. "dist". A n_vertices x n_neighbors matrix containing the distances of the nearest neighbors.
#' @param n_trees Number of trees to build when constructing the nearest neighbor index. The more trees specified, the larger the index, but the better the results. With search_k, determines the accuracy of the Annoy nearest neighbor search. Only used if the nn_method is "annoy". Sensible values are between 10 to 100.
#' @param search_k Number of nodes to search during the neighbor retrieval. The larger k, the more the accurate results, but the longer the search takes. With n_trees, determines the accuracy of the Annoy nearest neighbor search. Only used if the nn_method is "annoy".
#' @param approx_pow If TRUE, use an approximation to the power function in the UMAP gradient.
#' @param y Optional target data for supervised dimension reduction. Can be a vector, matrix or data frame. Use the target_metric parameter to specify the metrics to use, using the same syntax as metric. Usually either a single numeric or factor column is used, but more complex formats are possible.
#' @param target_n_neighbors Number of nearest neighbors to use to construct the target simplicial set. Default value is n_neighbors. Applies only if y is non-NULL and numeric.
#' @param target_metric The metric used to measure distance for y if using supervised dimension reduction. Used only if y is numeric.
#' @param target_weight Weighting factor between data topology and target topology. A value of 0.0 weights entirely on data, a value of 1.0 weights entirely on target. The default of 0.5 balances the weighting equally between data and target. Only applies if y is non-NULL.
#' @param pca f set to a positive integer value, reduce data to this number of columns using PCA. Doesn't applied if the distance metric is "hamming", or the dimensions of the data is larger than the number specified (i.e. number of rows and columns must be larger than the value of this parameter). If you have > 100 columns in a data frame or matrix, reducing the number of columns in this way may substantially increase the performance of the nearest neighbor search at the cost of a potential decrease in accuracy. In many t-SNE applications, a value of 50 is recommended, although there's no guarantee that this is appropriate for all settings.
#' @param pca_center If TRUE, center the columns of X before carrying out PCA. For binary data, it's recommended to set this to FALSE.
#' @param pcg_rand If TRUE, use the PCG random number generator (O'Neill, 2014) during optimization. Otherwise, use the faster (but probably less statistically good) Tausworthe "taus88" generator. The default is TRUE.
#' @param fast_sgd If TRUE, then the following combination of parameters is set: pcg_rand = TRUE, n_sgd_threads = "auto" and approx_pow = TRUE. The default is FALSE. Setting this to TRUE will speed up the stochastic optimization phase, but give a potentially less accurate embedding, and which will not be exactly reproducible even with a fixed seed. For visualization, fast_sgd = TRUE will give perfectly good results. For more generic dimensionality reduction, it's safer to leave fast_sgd = FALSE. If fast_sgd = TRUE, then user-supplied values of pcg_rand, n_sgd_threads, and approx_pow are ignored.
#' @param ret_model If TRUE, then return extra data that can be used to add new data to an existing embedding via umap_transform. The embedded coordinates are returned as the list item embedding. If FALSE, just return the coordinates. This parameter can be used in conjunction with ret_nn and ret_extra. Note that some settings are incompatible with the production of a UMAP model: external neighbor data (passed via a list to nn_method), and factor columns that were included via the metric parameter. In the latter case, the model produced is based only on the numeric data. A transformation using new data is possible, but the factor columns in the new data are ignored.
#' @param ret_nn If TRUE, then in addition to the embedding, also return nearest neighbor data that can be used as input to nn_method to avoid the overhead of repeatedly calculating the nearest neighbors when manipulating unrelated parameters (e.g. min_dist, n_epochs, init). See the "Value" section for the names of the list items. If FALSE, just return the coordinates. Note that the nearest neighbors could be sensitive to data scaling, so be wary of reusing nearest neighbor data if modifying the scale parameter. This parameter can be used in conjunction with ret_model and ret_extra.
#' @param n_threads Number of threads to use (except during stochastic gradient descent). Default is half the number of concurrent threads supported by the system. For nearest neighbor search, only applies if nn_method = "annoy". If n_threads > 1, then the Annoy index will be temporarily written to disk in the location determined by tempfile.
#' @param n_sgd_threads Number of threads to use during stochastic gradient descent. If set to > 1, then be aware that if batch = FALSE, results will not be reproducible, even if set.seed is called with a fixed seed before running. Set to "auto" to use the same value as n_threads.
#' @param grain_size The minimum amount of work to do on each thread. If this value is set high enough, then less than n_threads or n_sgd_threads will be used for processing, which might give a performance improvement if the overhead of thread management and context switching was outweighing the improvement due to concurrent processing. This should be left at default (1) and work will be spread evenly over all the threads specified.
#' @param tmpdir Temporary directory to store nearest neighbor indexes during nearest neighbor search. Default is tempdir. The index is only written to disk if n_threads > 1 and nn_method = "annoy"; otherwise, this parameter is ignored.
#' @param verbose If TRUE, log details to the console.
#'
#' @return An object of class iCellR.
#' @import uwot
#' @export
run.umap <- function (x = NULL,
                      my.seed = 0,
                      dims = 1:10,
                      n_neighbors = 15,
                      n_components = 2,
                      metric = "euclidean",
                      n_epochs = NULL,
                      learning_rate = 1,
                      scale = FALSE,
                      init = "spectral",
                      init_sdev = NULL,
                      spread = 1,
                      min_dist = 0.01,
                      set_op_mix_ratio = 1,
                      local_connectivity = 1,
                      bandwidth = 1,
                      repulsion_strength = 1,
                      negative_sample_rate = 5,
                      a = NULL,
                      b = NULL,
                      nn_method = NULL,
                      n_trees = 50,
                      search_k = 2 * n_neighbors * n_trees,
                      approx_pow = FALSE,
                      y = NULL,
                      target_n_neighbors = n_neighbors,
                      target_metric = "euclidean",
                      target_weight = 0.5,
                      pca = NULL,
                      pca_center = TRUE,
                      pcg_rand = TRUE,
                      fast_sgd = FALSE,
                      ret_model = FALSE,
                      ret_nn = FALSE,
                      n_threads = 1,
                      n_sgd_threads = 0,
                      grain_size = 1,
                      tmpdir = tempdir(),
                      verbose = getOption("verbose", TRUE)) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  # https://github.com/lmcinnes/umap
  # get PCA data
  set.seed(my.seed)
  DATA <- x@pca.data
  DATA <- DATA[dims]
################ uwot
  myUMAP = umap(DATA, n_neighbors = n_neighbors,
                n_components = n_components,
                metric = metric,
                n_epochs = n_epochs,
                learning_rate = learning_rate,
                scale = scale,
                init = init,
                init_sdev = init_sdev,
                spread = spread,
                min_dist = min_dist,
                set_op_mix_ratio = set_op_mix_ratio,
                local_connectivity = local_connectivity,
                bandwidth = bandwidth,
                repulsion_strength = repulsion_strength,
                negative_sample_rate = negative_sample_rate,
                a = a,
                b = b,
                nn_method = nn_method,
                n_trees = n_trees,
                search_k = search_k,
                approx_pow = approx_pow,
                y = y,
                target_n_neighbors = target_n_neighbors,
                target_metric = target_metric,
                target_weight = target_weight,
                pca = pca,
                pca_center = pca_center,
                pcg_rand = pcg_rand,
                fast_sgd = fast_sgd,
                ret_model = ret_model,
                ret_nn = ret_nn,
                n_threads = n_threads,
                n_sgd_threads = n_sgd_threads,
                grain_size = grain_size,
                tmpdir = tmpdir,
                verbose = verbose)
#############
  myUMAP = as.data.frame((myUMAP))
  My.distances = as.data.frame(as.matrix(dist((DATA))))[1]
  colnames(My.distances) <- "V3"
  myUMAP <- cbind(myUMAP,My.distances)
  attributes(x)$umap.data <- myUMAP
# return
  return(x)
}
