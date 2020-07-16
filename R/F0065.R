#' Run anchor alignment on the main data.
#'
#' This function takes an object of class iCellR and runs anchor alignment. It's a wrapper for Seurat.
#' @param x An object of class iCellR.
#' @param method Choose from "base.mean.rank" or "gene.model", default is "base.mean.rank". If gene.model is chosen you need to provide gene.list.
#' @param top.rank A number taking the top genes ranked by base mean, default = 500.
#' @param data.type Choose from "main" and "imputed", default = "main"
#' @param gene.list A charactor vector of genes to be used for PCA. If "clust.method" is set to "gene.model", default = "my_model_genes.txt".
#' @param normalization.method Choose from "LogNormalize", "CLR" and "RC". LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p. CLR: Applies a centered log ratio transformation. RC: Relative counts. Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. No log-transformation is applied. For counts per million (CPM) set ‘scale.factor = 1e6’
#' @param scale.factor Sets the scale factor for cell-level normalization.
#' @param margin If performing CLR normalization, normalize across features (1) or cells (2)
#' @param block.size How many cells should be run in each chunk, will try to split evenly across threads
#' @param selection.method Choose from "vst","mean.var.plot (mvp)","dispersion (disp)".
#' @param nfeatures Number of features to select as top variable features; only used when ‘selection.method’ is set to ‘'dispersion'’ or ‘'vst'’
#' @param anchor.features A numeric value. This will call ‘SelectIntegrationFeatures’ to select the provided number of features to be used in anchor finding
#' @param scale Whether or not to scale the features provided. Only set to FALSE if you have previously scaled the features you want to use for each object in the object.list
#' @param sct.clip.range Numeric of length two specifying the min and max values the Pearson residual will be clipped to
#' @param reduction cca: Canonical correlation analysis. rpca: Reciprocal PCA
#' @param l2.norm Perform L2 normalization on the CCA cell embeddings after dimensional reduction
#' @param dims Which dimensions to use from the CCA to specify the neighbor search space
#' @param k.anchor How many neighbors (k) to use when picking anchors
#' @param k.filter How many neighbors (k) to use when filtering anchors
#' @param k.score How many neighbors (k) to use when scoring anchors
#' @param max.features The maximum number of features to use when specifying the neighborhood search space in the anchor filtering
#' @param nn.method Method for nearest neighbor finding. Options include: rann, annoy
#' @param eps Error bound on the neighbor finding algorithm (from RANN)
#' @param k.weight Number of neighbors to consider when weighting
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
                     margin = 1,
                     block.size = NULL,
                     selection.method = "vst",
                     nfeatures = 2000,
                     anchor.features = 2000,
                     scale = TRUE,
                     sct.clip.range = NULL,
                     reduction = c("cca", "rpca"),
                     l2.norm = TRUE,
                     dims = 1:30,
                     k.anchor = 5,
                     k.filter = 200,
                     k.score = 30,
                     max.features = 200,
                     nn.method = "rann",
                     eps = 0,
                     k.weight = 100) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  #
  ##########
  start_time1 <- Sys.time()
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
    IDs = grep(paste("^",i,sep=""), Cells, value = TRUE)
    mydata <- DATA[ , which(names(DATA) %in% IDs)]
    MyMassg <- paste("    Preparing sample:",as.character(as.matrix(strsplit(i,'_',fixed=TRUE))))
    message(MyMassg)
    mydata <- CreateSeuratObject(counts = mydata, project = i)
    mydata <- NormalizeData(mydata, normalization.method = "LogNormalize",
                            scale.factor =scale.factor,
                            margin = margin, block.size=block.size)
    mydata <- FindVariableFeatures(mydata,
                                   selection.method = selection.method,
                                   nfeatures = nfeatures)
    SampNam <- paste("iCellRSample",i,sep="_")
    eval(call("<-", as.name(SampNam), mydata))
  }
  ######## get objects
  filenames <- ls(pattern="iCellRSample_")
  object.list <- mget(filenames)
#####
  if (normalization.method == "SCT") {
    for (i in names(object.list)) {
      object.list[[i]] <- SCTransform(object.list[[i]], verbose = TRUE)
    }
    anchor.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = nfeatures)
    object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = anchor.features)
  }
######
  Myanchors <- FindIntegrationAnchors(object.list = object.list,
                                      anchor.features = anchor.features, scale = scale,
                                      normalization.method = normalization.method,
                                      sct.clip.range = sct.clip.range, reduction = reduction,
                                      l2.norm = l2.norm, dims = dims, k.anchor = k.anchor, k.filter = k.filter,
                                      k.score = k.score, max.features = max.features, nn.method = nn.method, eps = eps)
  MYcombined <- IntegrateData(anchorset = Myanchors, dims = dims,
                              normalization.method = normalization.method,
                              k.weight = k.weight)
  ####
  DefaultAssay(MYcombined) <- "integrated"
  MYcombined <- ScaleData(MYcombined)
#  MYcombined <- RunPCA(MYcombined, npcs = 30)
  data <- as.data.frame(MYcombined@assays$integrated@scale.data)
  ############
  BestOrder <- colnames(DATA)
  #### order colnames by the original data
  data <- data[ , order(match(names(data),BestOrder))]
  ### PCA
  ###
  message(" Running PCA ...")
  counts.pca <- prcomp(data, center = TRUE, scale. = FALSE)
  attributes(x)$pca.info <- counts.pca
  dataPCA = data.frame(counts.pca$rotation) # [1:max.dim]
  attributes(x)$pca.data <- dataPCA
  #########
  end_time1 <- Sys.time()
  Time = difftime(end_time1,start_time1,units = "mins")
  Time = round(as.numeric(Time),digits = 2)
  message(paste("Total time",Time,"mins"))
  message("All done!")
  return(x)
}
