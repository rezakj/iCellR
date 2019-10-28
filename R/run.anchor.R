#' Run anchor alignment on the main data.
#'
#' This function takes an object of class iCellR and runs anchor alignment. It's a wrapper for Seurat.
#' @param x An object of class iCellR.
#' @param method Choose from "base.mean.rank" or "gene.model", default is "base.mean.rank". If gene.model is chosen you need to provide gene.list.
#' @param top.rank A number taking the top genes ranked by base mean, default = 500.
#' @param data.type Choose from "main" and "imputed", default = "main"
#' @param gene.list A charactor vector of genes to be used for PCA. If "clust.method" is set to "gene.model", default = "my_model_genes.txt".
#' @return An object of class iCellR.
#' @importFrom htmlwidgets saveWidget
#' @export
run.anchor <- function (x = NULL,
                     method = "base.mean.rank",
                     top.rank = 500,
                     gene.list = "character",
                     data.type = "main",
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     margin = 1, block.size = NULL,
                     selection.method = "vst",
                     nfeatures = 2000,
                     anchor.features = 2000,
                     scale = TRUE,
                     sct.clip.range = NULL, reduction = c("cca", "rpca"),
                     l2.norm = TRUE, dims = 1:30, k.anchor = 5, k.filter = 200,
                     k.score = 30, max.features = 200, nn.method = "rann", eps = 0,
                     k.weight = 100,
                     verbose = TRUE) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  #
  ##########
  #  require(Seurat)
  # Get data
  ## get main data
  if (data.type == "main") {
    DATA <- x@main.data
  }
  if (data.type == "imputed") {
    DATA <- x@imputed.data
  }
  ####
  if (method == "base.mean.rank") {
    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = TRUE), ]
    DATA <- head(raw.data.order,top.rank)
  }
  # gene model
  if (method == "gene.model") {
    if (gene.list[1] == "character") {
      stop("please provide gene names for clustering")
    } else {
      genesForClustering <- gene.list
      DATA <- subset(DATA, rownames(DATA) %in% genesForClustering)
    }
  }
  # Get conds
  Cells <- colnames(x@main.data)
  Conds <- as.character((unique(data.frame(do.call('rbind', strsplit(as.character(Cells),'_',fixed=TRUE)))[1]))$X1)
  ############
  if (length(Conds) == 1) {
    stop("You need more then one condition/sample to run this function")
  }
  ## get data
  Patt <- paste(Conds, "_",sep="")
  ################
  ###########
  if(!"Seurat" %in% (.packages())){
    stop("Please load Seurat package: library(Seurat)")
  }
  ##########
  message(" Preparing samples ...")
  for(i in Patt){
    IDs = grep(i, Cells, value = TRUE)
    mydata <- DATA[ , which(names(DATA) %in% IDs)]
    MyMassg <- paste("    Preparing sample:",as.character(as.matrix(strsplit(i,'_',fixed=TRUE))))
    message(MyMassg)
    mydata <- CreateSeuratObject(counts = mydata, project = i)
    mydata <- NormalizeData(mydata, normalization.method=normalization.method,
                            scale.factor =scale.factor, verbose = verbose,
                            margin = margin, block.size=block.size)
    mydata <- FindVariableFeatures(mydata,
                                   selection.method = selection.method,
                                   nfeatures = nfeatures)
    SampNam <- paste("iCellRSample",i,sep="_")
    eval(call("<-", as.name(SampNam), mydata))
  }
  ######## get objects
  filenames <- ls(pattern="iCellRSample_")
  Myanchors <- FindIntegrationAnchors(object.list = mget(filenames),
                                      anchor.features = anchor.features, scale = scale,
                                      normalization.method = normalization.method,
                                      sct.clip.range = sct.clip.range, reduction = reduction,
                                      l2.norm = l2.norm, dims = dims, k.anchor = k.anchor, k.filter = k.filter,
                                      k.score = k.score, max.features = max.features, nn.method = nn.method, eps = eps)
  MYcombined <- IntegrateData(anchorset = Myanchors, dims = dims,
                              normalization.method = normalization.method,
                              verbose =verbose, k.weight = k.weight)
  DefaultAssay(MYcombined) <- "integrated"
  MYcombined <- ScaleData(MYcombined, verbose = FALSE)
#  MYcombined <- RunPCA(MYcombined, npcs = 30, verbose = FALSE)
  data <- as.data.frame(MYcombined@assays$integrated@scale.data)
  ############
  BestOrder <- colnames(DATA)
  #### order colnames by the original data
  data <- data[ , order(match(names(data),BestOrder))]
  ### PCA
  message(" Running PCA ...")
  counts.pca <- prcomp(data, center = TRUE, scale. = FALSE)
  attributes(x)$pca.info <- counts.pca
  dataPCA = data.frame(counts.pca$rotation) # [1:max.dim]
  message("All done!")
  attributes(x)$pca.data <- dataPCA
  return(x)
}
