#' Run PCA on the main data
#'
#' This function takes an object of class iCellR and runs PCA on the main data.
#' @param x An object of class iCellR.
#' @param clust.method Choose from "base.mean.rank" or "gene.model", default is "base.mean.rank".
#' @param top.rank A number taking the top genes ranked by base mean, default = 500.
#' @param plus.log.value A number to add to each value in the matrix before log transformasion to aviond Inf numbers, default = 0.1.
#' @param gene.list A list of genes to be used for PCA. If "clust.method" is set to "gene.model", default = "my_model_genes.txt".
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- run.pca(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")
#' }
#' @export
run.cca <- function (x = NULL,
                     top.vari.genes = 1000,
                     cc.number = 30,
                     dims.align = 1:20,
                     normalize.data = TRUE,
                     scale.data = TRUE,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     display.progress = TRUE) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  #
  require(Seurat)
  # Get data
  DATA <- x@main.data
  # Get conds
  Cells <- colnames(x@main.data)
  Conds <- as.character((unique(data.frame(do.call('rbind', strsplit(as.character(Cells),'_',fixed=TRUE)))[1]))$X1)
  if (length(Conds) == 1) {
    stop("You need more then one condition/sample to run this function")
  }
## get data
    Patt <- paste(Conds[1], "_",sep="")
    FistCond = grep(Patt, Cells, value = T)
    DATA1 <- DATA[ , which(names(DATA) %in% FistCond)]
    DATA2 <- DATA[ , -which(names(DATA) %in% FistCond)]
    # make Seurat objects
    # Set up DATA1
    ctrl <- CreateSeuratObject(raw.data = DATA1, project = "IMMUNE_CTRL")
    ctrl@meta.data$stim <- "CTRL"
    if (normalize.data == TRUE) {
      ctrl <- NormalizeData(ctrl, normalization.method = normalization.method,
                            scale.factor = scale.factor,
                            display.progress = display.progress)
    }
    if (scale.data == TRUE) {
      ctrl <- ScaleData(ctrl, display.progress = display.progress)
    }
    ctrl <- FindVariableGenes(ctrl, do.plot = F)
    # Set up DATA2
    stim <- CreateSeuratObject(raw.data = DATA2, project = "IMMUNE_STIM")
    stim@meta.data$stim <- "STIM"
    if (normalize.data == TRUE) {
      stim <- NormalizeData(stim, normalization.method = normalization.method,
                            scale.factor = scale.factor,
                            display.progress = display.progress)
    }
    if (scale.data == TRUE) {
      stim <- ScaleData(stim, display.progress = display.progress)
    }
    stim <- FindVariableGenes(stim, do.plot = F)
    # get variable genes
    g.1 <- head(rownames(ctrl@hvg.info), top.vari.genes)
    g.2 <- head(rownames(stim@hvg.info), top.vari.genes)
    genes.use <- unique(c(g.1, g.2))
    genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
    genes.use <- intersect(genes.use, rownames(stim@scale.data))
    # Runn CCA
    combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = cc.number)
    # alignement
    combined <- AlignSubspace(combined,
                                     reduction.type = "cca",
                                     grouping.var = "stim",
                                     dims.align = dims.align)
    # Get CCA
    ToCCA <- as.data.frame(combined@dr$cca.aligned@cell.embeddings)
  # object
  attributes(x)$cca.data <- ToCCA
  return(x)
}
