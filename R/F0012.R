#' Run PCA on the main data
#'
#' This function takes an object of class iCellR and runs PCA on the main data.
#' @param x An object of class iCellR.
#' @param method Choose from "base.mean.rank" or "gene.model", default is "base.mean.rank". If gene.model is chosen you need to provide gene.list.
#' @param top.rank A number. Taking the top genes ranked by base mean, default = 500.
#' @param data.type Choose from "main" and "imputed", default = "main"
#' @param plus.log.value A number to add to each value in the matrix before log transformasion to aviond Inf numbers, default = 0.1.
#' @param gene.list A charactor vector of genes to be used for PCA. If "clust.method" is set to "gene.model", default = "my_model_genes.txt".
#' @param scale.data If TRUE the data will be scaled (log2 + plus.log.value), default = TRUE.
#' @return An object of class iCellR.
#' @examples
#' demo.obj <- run.pca(demo.obj, method = "gene.model", gene.list = demo.obj@gene.model)
#'
#' head(demo.obj@pca.data)[1:5]
#'
#' @export
run.pca <- function (x = NULL,
                          data.type = "main",
                          method = "base.mean.rank",
                          top.rank = 500,
                          plus.log.value = 0.1,
                          scale.data = TRUE,
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
    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = TRUE), ]
    TopNormLogScale <- head(raw.data.order,top.rank)
  }
  # gene model
  if (method == "gene.model") {
    if (gene.list[1] == "character") {
      stop("please provide gene names for clustering")
    } else {
      genesForClustering <- gene.list
      TopNormLogScale <- subset(DATA, rownames(DATA) %in% genesForClustering)
    }
  }
# Scale
  if(scale.data == TRUE) {
            TopNormLogScale <- log(TopNormLogScale + plus.log.value)
   }
# Returns
  # info
    counts.pca <- prcomp(TopNormLogScale, center = FALSE, scale. = FALSE)
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
