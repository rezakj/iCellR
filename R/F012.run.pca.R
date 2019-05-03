#' Run PCA on the main data
#'
#' This function takes an object of class iCellR and runs PCA on the main data.
#' @param x An object of class iCellR.
#' @param method Choose from "base.mean.rank" or "gene.model", default is "base.mean.rank". If gene.model is chosen you need to provide gene.list.
#' @param top.rank A number taking the top genes ranked by base mean, default = 500.
#' @param plus.log.value A number to add to each value in the matrix before log transformasion to aviond Inf numbers, default = 0.1.
#' @param gene.list A charactor vector of genes to be used for PCA. If "clust.method" is set to "gene.model", default = "my_model_genes.txt".
#' @param batch.norm If TRUE the data will be normalized based on the genes in gene.list or top ranked genes.
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- run.pca(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")
#' }
#' @export
run.pca <- function (x = NULL,
                          data.type = "main",
                          method = "base.mean.rank",
                          top.rank = 500,
                          plus.log.value = 0.1,
                          batch.norm = F,
                          gene.list = "character") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  # geth the genes and scale them based on model
  ## get main data
  if (data.type == "main") {
    DATA <- x@main.data
  }
  if (data.type == "imputed") {
    DATA <- x@imputed.data
  }
  # model base mean rank
  if (method == "base.mean.rank") {
    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = T), ]
    topGenes <- head(raw.data.order,top.rank)
    TopNormLogScale <- log(topGenes + plus.log.value)
    # TopNormLogScale <- scale(topGenes)
#    TopNormLogScale <- t(TopNormLogScale)
#    TopNormLogScale <- as.data.frame(t(scale(TopNormLogScale)))
  }
  # gene model
  if (method == "gene.model") {
    if (gene.list[1] == "character") {
      stop("please provide gene names for clustering")
    } else {
      genesForClustering <- gene.list
      topGenes <- subset(DATA, rownames(DATA) %in% genesForClustering)
      if (batch.norm == F){
         TopNormLogScale <- log(topGenes + plus.log.value)
        # TopNormLogScale <- scale(topGenes)
      }
      if (batch.norm == T){
        ## new method
        libSiz <- colSums(topGenes)
        norm.facts <- as.numeric(libSiz) / mean(as.numeric(libSiz))
        dataMat <- as.matrix(topGenes)
        normalized <- as.data.frame(sweep(dataMat, 2, norm.facts, `/`))
         TopNormLogScale <- log(normalized + plus.log.value)
        TopNormLogScale <- normalized
      }
    }
  }
# Returns
  # info
    counts.pca <- prcomp(TopNormLogScale, center = T, scale. = T)
    attributes(x)$pca.info <- counts.pca
    # DATA
    dataPCA = data.frame(counts.pca$rotation) # [1:max.dim]
    attributes(x)$pca.data <- dataPCA
    # optimal
    DATA <- counts.pca$sdev
    OPTpcs <- mean(DATA)*2
    OPTpcs <- (DATA > OPTpcs)
    OPTpcs <- length(OPTpcs[OPTpcs==TRUE]) + 1
    attributes(x)$opt.pcs <- OPTpcs
    # object
  return(x)
}
