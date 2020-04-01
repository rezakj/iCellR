#' Run MNN alignment on the main data.
#'
#' This function takes an object of class iCellR and runs MNN alignment. It's a wrapper for scran.
#' @param x An object of class iCellR.
#' @param method Choose from "base.mean.rank" or "gene.model", default is "base.mean.rank". If gene.model is chosen you need to provide gene.list.
#' @param top.rank A number taking the top genes ranked by base mean, default = 500.
#' @param data.type Choose from "main" and "imputed", default = "main"
#' @param gene.list A charactor vector of genes to be used for PCA. If "clust.method" is set to "gene.model", default = "my_model_genes.txt".
#' @param k An integer scalar specifying the number of nearest neighbors to consider when identifying MNNs.
#' @param cos.norm A logical scalar indicating whether cosine normalization should be performed on the input data prior to calculating distances between cells.
#' @param ndist A numeric scalar specifying the threshold beyond which neighbours are to be ignored when computing correction vectors. Each threshold is defined in terms of the number of median distances.
#' @param d Number of dimentions to pass to ‘multiBatchPCA’.
#' @param approximate Further arguments to pass to ‘multiBatchPCA’. Setting ‘approximate=TRUE’ is recommended for large data sets with many cells.
#' @param irlba.args Further arguments to pass to ‘multiBatchPCA’. Setting ‘approximate=TRUE’ is recommended for large data sets with many cells.
#' @param subset.row See ‘?"scran-gene-selection"’.
#' @param auto.order Logical scalar indicating whether re-ordering of batches should be performed to maximize the number of MNN pairs at each step. Alternatively an integer vector containing a permutation of ‘1:N’ where ‘N’ is the number of batches.
#' @param compute.variances Logical scalar indicating whether the percentage of variance lost due to non-orthogonality should be computed.
#' @param pc.input Logical scalar indicating whether the values in ‘...’ are already low-dimensional, e.g., the output of ‘multiBatchPCA’.
#' @param assay.type A string or integer scalar specifying the assay containing the expression values, if SingleCellExperiment objects are present in ‘...’.
#' @param get.spikes See ‘?"scran-gene-selection"’. Only relevant if ‘...’ contains SingleCellExperiment objects.
#' @param BNPARAM A BiocNeighborParam object specifying the nearest neighbor algorithm. Defaults to an exact algorithm if ‘NULL’, see ‘?findKNN’ for more details.
#' @param BPPARAM A BiocParallelParam object specifying whether the PCA and nearest-neighbor searches should be parallelized.
#' @return An object of class iCellR.
#' @importFrom htmlwidgets saveWidget
#' @export
run.mnn <- function (x = NULL,
                     method = "base.mean.rank",
                     top.rank = 500,
                     gene.list = "character",
                     data.type = "main",
                     k=20, cos.norm=TRUE, ndist=3, d=50, approximate=FALSE,
                     irlba.args=list(), subset.row=NULL, auto.order=FALSE, pc.input=FALSE,
                     compute.variances=FALSE, assay.type="logcounts", get.spikes=FALSE,
                     BNPARAM=NULL, BPPARAM=SerialParam()) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  #
  ##########
  start_time1 <- Sys.time()
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
  if(!"scran" %in% (.packages())){
    stop("Please load scran package: library(scran)")
  }
  ##########
  message(" Preparing samples ...")
  for(i in Patt){
    IDs = grep(paste("^",i,sep=""), Cells, value = TRUE)
    mydata <- DATA[ , which(names(DATA) %in% IDs)]
    CellN <- dim(mydata)[2]
    SampNam <- paste("iCellRSample",i,sep="_")
    mydata <- SingleCellExperiment(list(counts=as.matrix(mydata)))
    MyMassg <- paste("    Preparing sample:",as.character(as.matrix(strsplit(i,'_',fixed=TRUE))))
    message(MyMassg)
    #mydata <- computeSumFactors(mydata)
    # not normalize
    MyList <- list(rep(1,CellN))
    names(MyList) <- "size_factor"
    mydata@int_colData@listData <- MyList
    ####
    mydata <- normalize(mydata)
    eval(call("<-", as.name(SampNam), mydata))
  }
######## get objects
  filenames <- ls(pattern="iCellRSample_")
  for(i in filenames){
    mydata <- colnames(get(i))
    SampNam <- paste("MyColNames",i,sep="_")
    eval(call("<-", as.name(SampNam), mydata))
  }
  CellIDs <- as.character(unlist(mget(ls(pattern="MyColNames_"))))
##########
  filenames <- paste(filenames,collapse=",")
  ########## PARAMs
  MyPARAM1 <- paste("k=k,cos.norm=cos.norm,ndist=ndist,")
  MyPARAM2 <- paste("d=d,approximate=approximate,irlba.args=irlba.args,")
  MyPARAM3 <- paste("subset.row=subset.row,auto.order=auto.order,")
  MyPARAM4 <- paste("pc.input=pc.input,compute.variances=compute.variances,")
  MyPARAM5 <- paste("assay.type=assay.type,get.spikes=get.spikes,BNPARAM=BNPARAM,BPPARAM=BPPARAM)")
  MyPARAM <- paste(MyPARAM1,MyPARAM2,MyPARAM3,MyPARAM4,MyPARAM5)
  ############
  ha=paste("out <- fastMNN(",filenames,",",MyPARAM )
  message("   Running MNN ...")
  eval(parse(text=ha))
  MNN <- as.data.frame(t(out$corrected))
  ######### Col Names
  colnames(MNN) <- CellIDs
  #### Col Names order
  BestOrder <- colnames(DATA)
#### order colnames by the original data
  MNN <- MNN[ , order(match(names(MNN),BestOrder))]
  ### PCA
  message(" Running PCA ...")
  counts.pca <- prcomp((MNN), center = TRUE, scale. = FALSE)
  attributes(x)$pca.info <- counts.pca
  dataPCA = data.frame(counts.pca$rotation) # [1:max.dim]
  attributes(x)$pca.data <- dataPCA
#########
  message("All done!")
  end_time1 <- Sys.time()
  Time = difftime(end_time1,start_time1,units = "mins")
  Time = round(as.numeric(Time),digits = 2)
  message(paste("Total time",Time,"mins"))
  return(x)
}
