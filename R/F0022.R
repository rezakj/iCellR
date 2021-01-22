#' Impute the main data
#'
#' This function takes an object of class iCellR and runs imputation on the main data.
#' @param x An object of class iCellR.
#' @param imp.method Choose between "iCellR.imp" and "magic", defualt = "iCellR.imp".
#' @param nn Number of neighboring cells to find, default = 10.
#' @param dims PC dimentions to be used for the analysis, default = 10.
#' @param data.type Choose between "tsne", "pca", "umap", "diffusion", "knetl", default = "pca".
#' @param genes character or integer vector, default: NULL vector of column names or column indices for which to return smoothed data If 'all_genes' or NULL, the entire smoothed matrix is returned
#' @param k if imp.method is magic; int, optional, default: 10 number of nearest neighbors on which to build kernel
#' @param alpha if imp.method is magic; int, optional, default: 15 sets decay rate of kernel tails. If NULL, alpha decaying kernel is not used
#' @param t if imp.method is magic; int, optional, default: 'auto' power to which the diffusion operator is powered sets the level of diffusion. If 'auto', t is selected according to the Procrustes disparity of the diffused data.'
#' @param npca number of PCA components that should be used; default: 100.
#' @param init magic object, optional object to use for initialization. Avoids recomputing intermediate steps if parameters are the same.
#' @param t.max if imp.method is magic; int, optional, default: 20 Maximum value of t to test for automatic t selection.
#' @param knn.dist.method string, optional, default: 'euclidean'. recommended values: 'euclidean', 'cosine' Any metric from 'scipy.spatial.distance' can be used distance metric for building kNN graph.
#' @param verbose 'int' or 'boolean', optional (default : 1) If 'TRUE' or '> 0', message verbose updates.
#' @param n.jobs 'int', optional (default: 1) The number of jobs to use for the computation. If -1 all CPUs are used. If 1 is given, no parallel computing code is used at all, which is useful for debugging. For n_jobs below -1, (n.cpus + 1 + n.jobs) are used. Thus for n_jobs = -2, all CPUs but one are used
#' @param seed int or 'NULL', random state (default: 'NULL')
#' @return An object of class iCellR.
#' @import progress
#' @export
run.impute <- function (x = NULL,
                        imp.method = "iCellR.imp", dims = 1:10, nn = 10,
                        data.type = "pca",genes = "all_genes", k = 10, alpha = 15, t = "auto",
                        npca = 100, init = NULL, t.max = 20,
                        knn.dist.method = "euclidean", verbose = 1, n.jobs = 1,
                        seed = NULL) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  # get data
  start_time1 <- Sys.time()
  #####
  DATA <- x@main.data
  message(paste(" main data dimentions:",dim(DATA)[1],"genes and",dim(DATA)[2],"cells"))
  message(" Genes with no coverage are being removed from the matrix ...")
  DATA <- DATA[ rowSums(DATA) > 0, ]
  message(paste(" dimentions after removing zero cov. genes:",dim(DATA)[1],"genes and",dim(DATA)[2],"cells"))
  #########
  #########################
  if (imp.method == "iCellR.imp"){
    my.data = DATA
    if(data.type == "pca") {
      my.data.my.pca= (t(x@pca.data))[dims, ]
#      num.of.PCs = c(dims)
#      my.data.my.pca = t(x@pca.info$rotation)[num.of.PCs, ]
    }
    if(data.type == "umap") {
      my.data.my.pca = t(x@umap.data)
    }
    if(data.type == "tsne") {
      my.data.my.pca = t(x@tsne.data)
    }
    if(data.type == "knetl") {
      my.data.my.pca = t(x@knetl.data)
    }
    if(data.type == "diffusion") {
      my.data.my.pca = t(x@diffusion.data)
    }
#########
    start_time <- Sys.time()
    message(paste("   Calculating distance ..."))
    My.distances = as.matrix(dist(t(my.data.my.pca), method = knn.dist.method))
    end_time <- Sys.time()
    Time = difftime(end_time,start_time,units = "mins")
    Time = round(as.numeric(Time),digits = 2)
    message(paste("   Calculated distance in",Time,"mins"))
    ncells = dim(my.data)[2]
#    cell.num = ceiling(cell.ratio/100 * ncells)
    cell.num = nn
#    message(paste("    ",cell.ratio,"percent of ",ncells, "cells is", cell.num))
    message(paste("    Finding",cell.num, "neighboring cells per cell ..."))
    message("     To change the number of neighboring cells change it using the nn option")
    ### time
    pb <- progress_bar$new(total = ncells,
                           format = "[:bar] :current/:total (:percent) :elapsedfull eta: :eta",
                           clear = FALSE, width= 60)
#######
    KNN1 = lapply(1:ncells, function(findKNN){
      pb$tick()
      order(My.distances[,findKNN])[1:cell.num]})
    ############
    message(paste("   correcting the coverage of the neighboring cells (mean) ..."))
    ### time
#    pb <- progress_bar$new(total = ncells,
#                           format = "[:bar] :current/:total (:percent) :elapsedfull eta: :eta",
#                           clear = FALSE, width= 60)
    ##########
#    data.sum1 = sapply(KNN1, function(sum.cov){
#      pb$tick()
#      rowMeans(my.data[, sum.cov])})
    #######
    my.data <- as.matrix(my.data)
    pb <- progress_bar$new(total = ncells,
                           format = "[:bar] :current/:total (:percent) :elapsedfull eta: :eta",
                           clear = FALSE, width= 60)
    ##########
    data.sum1 = sapply(KNN1, function(sum.cov){
      pb$tick()
      rowMeans(my.data[, sum.cov])})
#    GETmean(my.data[, sum.cov])})
    ############
    #my.data <- as.data.frame(my.data)
    data.sum1 <- as.data.frame(data.sum1)
    row.names(data.sum1) <- row.names(my.data)
    colnames(data.sum1) <- colnames(my.data)
    data.sum1 <- round(data.sum1, digits = 3)
######
    message(paste("All done!"))
    end_time1 <- Sys.time()
    Time = difftime(end_time1,start_time1,units = "mins")
    Time = round(as.numeric(Time),digits = 3)
    message(paste("Total time",Time,"mins"))
    attributes(x)$imputed.data <- data.sum1
    return(x)
  }
  ###############################################################
  if (imp.method == "magic"){
    ###########
    if(!"phateR" %in% (.packages())){
      stop("Please load phateR package: library(phateR). This function requires 'phateR','Rmagic'and'viridis'.")
    }
    ##########
    ###########
    if(!"Rmagic" %in% (.packages())){
      stop("Please load Rmagic package: library(Rmagic). This function requires 'phateR','Rmagic'and'viridis'.")
    }
    ##########
    ###########
    if(!"viridis" %in% (.packages())){
      stop("Please load viridis package: library(viridis). This function requires 'phateR','Rmagic'and'viridis'.")
    }
    ##########
    ### load packages
    # https://github.com/KrishnaswamyLab/MAGIC/tree/master/Rmagic
    # https://www.analyticsvidhya.com/blog/2016/03/tutorial-powerful-packages-imputing-missing-values/
    DATA <- as.data.frame(t(DATA))
    data_MAGIC <- magic(DATA,
                        genes = genes, k = k, alpha = alpha, t = t,
                        npca = npca, init = init, t.max = t.max,
                        knn.dist.method = knn.dist.method, verbose = verbose, n.jobs = n.jobs,
                        seed = seed)
    DATA <- as.data.frame(t(data_MAGIC$result))
    # return
    DATA <- round(DATA, digits = 3)
    attributes(x)$imputed.data <- DATA
    # return
    return(x)
  }
}
