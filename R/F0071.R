#' Cell cycle phase prediction
#'
#' This function takes an object of class iCellR and assignes cell cycle stage for the cells.
#' @param object A data frame containing gene counts for cells.
#' @param data.type Choose from "raw.data" or "main.data", "imputed.data", default = "main.data".
#' @param scoring.List Genes that are used as a marker for phases.
#' @param return.stats Return the data or object. If FALSE the object would be returned.
#' @param scoring.method Choose from "tirosh (Tirosh, et. al. 2016), mean, sum, gsva, ssgsea, zscore and plage. , default = "tirosh".
#' @return The data frame object
#' @importFrom Hmisc cut2
#' @export
i.score <- function (object = NULL,
                        data.type = "main.data",
                        scoring.List = NULL,
                        return.stats = TRUE,
                        scoring.method = "tirosh") {
  if ("iCellR" != class(object)[1]) {
    stop("object should be an object of class iCellR")
  }
  if (is.null(scoring.List)) {
    stop('At least 2 set of genes are required for scoring.List
         Example: scoring.List = c("G0","G1S","G2M","M","MG1","S")')
  }
  ##### get data
  ## get main data
  DATAraw <- object@raw.data
  ######
  if (data.type == "main.data") {
    DATA <- object@main.data
    all.genes <- row.names(DATA)
    object@raw.data <- object@main.data
  }
  if (data.type == "raw.data") {
    DATA <- object@raw.data
    all.genes <- row.names(DATA)
  }
  if (data.type == "imputed.data") {
    DATA <- object@imputed.data
    all.genes <- row.names(DATA)
    object@raw.data <- object@imputed.data
  }
  ####
  if (scoring.method == "tirosh") {
    object <- qc.stats(object)
    object <- cc(object)
  }
   #########
  if (scoring.method == "sum") {
    for (i in scoring.List) {
      data <- get(i)
      message(i)
      message(paste(" Finding",length(data),i,"genes ..."))
      data <- paste("^",data,"$", sep="")
      data <- paste(data,collapse="|")
      data <- grep(data, x = all.genes, value = TRUE, ignore.case = TRUE)
      message(paste("         Found",length(data),"genes for",i,"."))
      message("")
      subDATA <- subset(DATA,rownames(DATA) %in% data)
#      dim(subDATA)
#      subDATA <- as.data.frame(colMeans(subDATA))
      subDATA <- as.data.frame(colSums(subDATA))
      colnames(subDATA) <- i
      ###
      NameCol=paste("cellCyclePhase",i,sep="_")
      eval(call("<-", as.name(NameCol), subDATA))
    }
    ###### make data frame
    filenames <- ls(pattern="cellCyclePhase_")
    datalist <- mget(filenames)
    FinalData1 <- as.data.frame(t(as.data.frame(datalist)))
    ###
    FinalData <- hto.anno(hto.data = FinalData1,
                          cov.thr = 10,
                          assignment.thr = 50)
  }
  ########
  if (scoring.method == "mean") {
    for (i in scoring.List) {
      data <- get(i)
      message(i)
      message(paste(" Finding",length(data),i,"genes ..."))
      data <- paste("^",data,"$", sep="")
      data <- paste(data,collapse="|")
      data <- grep(data, x = all.genes, value = TRUE, ignore.case = TRUE)
      message(paste("         Found",length(data),"genes for",i,"."))
      message("")
      subDATA <- subset(DATA,rownames(DATA) %in% data)
      #      dim(subDATA)
      subDATA <- as.data.frame(colMeans(subDATA))
#      subDATA <- as.data.frame(colSums(subDATA))
      colnames(subDATA) <- i
      ###
      NameCol=paste("cellCyclePhase",i,sep="_")
      eval(call("<-", as.name(NameCol), subDATA))
    }
    ###### make data frame
    filenames <- ls(pattern="cellCyclePhase_")
    datalist <- mget(filenames)
    FinalData1 <- as.data.frame(t(as.data.frame(datalist)))
    ###
    FinalData <- hto.anno(hto.data = FinalData1,
                          cov.thr = 10,
                          assignment.thr = 50)
  }
  ########
  if (scoring.method == "ssgsea") {
    if(!"GSVA" %in% (.packages())){
      stop("Please load scran package: library(GSVA)")
    }
    for (i in scoring.List) {
      data <- get(i)
      message(i)
      message(paste(" Finding",length(data),i,"genes ..."))
      data <- paste("^",data,"$", sep="")
      data <- paste(data,collapse="|")
      data <- grep(data, x = all.genes, value = TRUE, ignore.case = TRUE)
      message(paste("         Found",length(data),"genes for",i,"."))
      message("")
      subDATA <- subset(DATA,rownames(DATA) %in% data)
      Signatures <- list(rownames(subDATA))
      Table <- as.matrix(DATA)
      Scores <- as.data.frame(t(gsva(Table, Signatures,method="ssgsea", kcdf="none", mx.diff=TRUE)))
      colnames(Scores) <- i
      ###
      NameCol=paste("cellCyclePhase",i,sep="_")
      eval(call("<-", as.name(NameCol), Scores))
    }
    ###### make data frame
    filenames <- ls(pattern="cellCyclePhase_")
    datalist <- mget(filenames)
    FinalData1 <- as.data.frame(t(as.data.frame(datalist)))
    ###
    FinalData <- hto.anno(hto.data = FinalData1,
                          cov.thr = 10,
                          assignment.thr = 50)
  }
  ########
  ########
  if (scoring.method == "gsva") {
    if(!"GSVA" %in% (.packages())){
      stop("Please load scran package: library(GSVA)")
    }
    for (i in scoring.List) {
      data <- get(i)
      message(i)
      message(paste(" Finding",length(data),i,"genes ..."))
      data <- paste("^",data,"$", sep="")
      data <- paste(data,collapse="|")
      data <- grep(data, x = all.genes, value = TRUE, ignore.case = TRUE)
      message(paste("         Found",length(data),"genes for",i,"."))
      message("")
      subDATA <- subset(DATA,rownames(DATA) %in% data)
      Signatures <- list(rownames(subDATA))
      Table <- as.matrix(DATA)
      Scores <- as.data.frame(t(gsva(Table, Signatures,method="gsva", kcdf="none", mx.diff=TRUE)))
      colnames(Scores) <- i
      ###
      NameCol=paste("cellCyclePhase",i,sep="_")
      eval(call("<-", as.name(NameCol), Scores))
    }
    ###### make data frame
    filenames <- ls(pattern="cellCyclePhase_")
    datalist <- mget(filenames)
    FinalData1 <- as.data.frame(t(as.data.frame(datalist)))
    ###
    FinalData <- hto.anno(hto.data = FinalData1,
                          cov.thr = 10,
                          assignment.thr = 50)
  }
  ########
  if (scoring.method == "zscore") {
    if(!"GSVA" %in% (.packages())){
      stop("Please load scran package: library(GSVA)")
    }
    for (i in scoring.List) {
      data <- get(i)
      message(i)
      message(paste(" Finding",length(data),i,"genes ..."))
      data <- paste("^",data,"$", sep="")
      data <- paste(data,collapse="|")
      data <- grep(data, x = all.genes, value = TRUE, ignore.case = TRUE)
      message(paste("         Found",length(data),"genes for",i,"."))
      message("")
      subDATA <- subset(DATA,rownames(DATA) %in% data)
      Signatures <- list(rownames(subDATA))
      Table <- as.matrix(DATA)
      Scores <- as.data.frame(t(gsva(Table, Signatures,method="zscore", kcdf="Gaussian", mx.diff=TRUE)))
      colnames(Scores) <- i
      ###
      NameCol=paste("cellCyclePhase",i,sep="_")
      eval(call("<-", as.name(NameCol), Scores))
    }
    ###### make data frame
    filenames <- ls(pattern="cellCyclePhase_")
    datalist <- mget(filenames)
    FinalData1 <- as.data.frame(t(as.data.frame(datalist)))
    ###
    FinalData <- hto.anno(hto.data = FinalData1,
                          cov.thr = 10,
                          assignment.thr = 50)
  }
  ########
  if (scoring.method == "plage") {
    if(!"GSVA" %in% (.packages())){
      stop("Please load scran package: library(GSVA)")
    }
    for (i in scoring.List) {
      data <- get(i)
      message(i)
      message(paste(" Finding",length(data),i,"genes ..."))
      data <- paste("^",data,"$", sep="")
      data <- paste(data,collapse="|")
      data <- grep(data, x = all.genes, value = TRUE, ignore.case = TRUE)
      message(paste("         Found",length(data),"genes for",i,"."))
      message("")
      subDATA <- subset(DATA,rownames(DATA) %in% data)
      Signatures <- list(rownames(subDATA))
      Table <- as.matrix(DATA)
      Scores <- as.data.frame(t(gsva(Table, Signatures,method="plage", kcdf="Gaussian", mx.diff=TRUE)))
      colnames(Scores) <- i
      ###
      NameCol=paste("cellCyclePhase",i,sep="_")
      eval(call("<-", as.name(NameCol), Scores))
    }
    ###### make data frame
    filenames <- ls(pattern="cellCyclePhase_")
    datalist <- mget(filenames)
    FinalData1 <- as.data.frame(t(as.data.frame(datalist)))
    ###
    FinalData <- hto.anno(hto.data = FinalData1,
                          cov.thr = 10,
                          assignment.thr = 50)
  }
  ########
  if (scoring.method == "tirosh") {
    ###############
    ##############
    ##############
    AddModuleScoreme <- function (object, genes.list = NULL, genes.pool = NULL, n.bin = 25,
                                  seed.use = 1, ctrl.size = 100, use.k = FALSE, enrich.name = "Cluster",
                                  random.seed = 1)
    {
      set.seed(seed = random.seed)
      genes.old <- genes.list
      if (use.k) {
        genes.list <- list()
        for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
          genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster ==
                                               i))
        }
        cluster.length <- length(x = genes.list)
      }
      else {
        if (is.null(x = genes.list)) {
          stop("Missing input gene list")
        }
        genes.list <- lapply(X = genes.list, FUN = function(x) {
          return(intersect(x = x, y = rownames(x = object@raw.data)))
        })
        cluster.length <- length(x = genes.list)
      }

      if (!all(LengthCheck(values = genes.list))) {
        warning(paste("Could not find enough genes in the object from the following gene lists:",
                      paste(names(x = which(x = !LengthCheck(values = genes.list)))),
                      "Attempting to match case..."))
        genes.list <- lapply(X = genes.old, FUN = CaseMatch,
                             match = rownames(x = object@raw.data))
      }
      if (!all(LengthCheck(values = genes.list))) {
        stop(paste("The following gene lists do not have enough genes present in the object:",
                   paste(names(x = which(x = !LengthCheck(values = genes.list)))),
                   "exiting..."))
      }

      if (is.null(x = genes.pool)) {
        genes.pool = rownames(x = object@raw.data)
      }
      data.avg <- Matrix::rowMeans(x = object@raw.data[genes.pool,
      ])
      data.avg <- data.avg[order(data.avg)]
      #
      data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg)/n.bin)))
      names(x = data.cut) <- names(x = data.avg)
      ctrl.use <- vector(mode = "list", length = cluster.length)
      for (i in 1:cluster.length) {
        genes.use <- genes.list[[i]]
        for (j in 1:length(x = genes.use)) {
          ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut ==
                                                                                  data.cut[genes.use[j]])], size = ctrl.size, replace = FALSE)))
        }
      }
      ctrl.use <- lapply(X = ctrl.use, FUN = unique)
      ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use),
                            ncol = ncol(x = object@raw.data))
      for (i in 1:length(ctrl.use)) {
        genes.use <- ctrl.use[[i]]
        ctrl.scores[i, ] <- Matrix::colMeans(x = object@raw.data[genes.use,
        ])
      }
      genes.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length,
                             ncol = ncol(x = object@raw.data))
      for (i in 1:cluster.length) {
        genes.use <- genes.list[[i]]
        data.use <- object@raw.data[genes.use, , drop = FALSE]
        genes.scores[i, ] <- Matrix::colMeans(x = data.use)
      }
      genes.scores.use <- genes.scores - ctrl.scores
      rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
      genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))
      rownames(x = genes.scores.use) <- colnames(x = object@raw.data)
      object <- AddMetaDatame(object = object, metadata = genes.scores.use,
                              col.name = colnames(x = genes.scores.use))
      gc(verbose = FALSE)
      return(object)
    }
    ########## 1
    LengthCheck <- function(values, cutoff = 0) {
      return(vapply(
        X = values,
        FUN = function(x) {
          return(length(x = x) > cutoff)
        },
        FUN.VALUE = logical(1)
      ))
    }
    #######
    ########## 2
    ####
    AddMetaDatame <- function (object, metadata, col.name = NULL)
    {
      if (typeof(x = metadata) != "list") {
        metadata <- as.data.frame(x = metadata)
        if (is.null(x = col.name)) {
          stop("Please provide a name for provided metadata")
        }
        colnames(x = metadata) <- col.name
      }
      cols.add <- colnames(x = metadata)
      meta.order <- match(rownames(object@stats), rownames(metadata))
      meta.add <- metadata[meta.order, ]
      if (all(is.null(x = meta.add))) {
        stop("Metadata provided doesn't match the cells in this object")
      }
      object@stats[, cols.add] <- meta.add
      return(object)
    }
    ##############
    ##############
    ##############
    for (i in scoring.List) {
      data <- get(i)
      message(i)
      message(paste(" Finding",length(data),i,"genes ..."))
      data <- paste("^",data,"$", sep="")
      data <- paste(data,collapse="|")
      data <- grep(data, x = all.genes, value = TRUE, ignore.case = TRUE)
      message(paste("         Found",length(data),"genes for",i,"."))
      message("")
      ###
      NameCol=paste("cellCyclePhase",i,sep="_")
      eval(call("<-", as.name(NameCol), data))
    }
    ########
    filenames <- ls(pattern="cellCyclePhase_")
    datalist <- mget(filenames)
    names(datalist) <- gsub("cellCyclePhase_","",names(datalist))
    #######
    enrich.name <- "Cell Cycle"
    genes.list <- datalist
    object.cc <- AddModuleScoreme(object = object, genes.list = genes.list,
                                  enrich.name = enrich.name,
                                  ctrl.size = min(vapply(X = genes.list,
                                                         FUN = length, FUN.VALUE = numeric(1))))
    #######
    TO <- length(names(datalist))+9
    FinalData1 <- (object.cc@stats)[10:TO]
    colnames(FinalData1) <- names(datalist)
    FinalData1 <- as.data.frame(t(FinalData1))
    ###
    FinalData <- hto.anno(hto.data = FinalData1,
                          cov.thr = 10,
                          assignment.thr = 50)
  }
  ########################
  if (return.stats == FALSE) {
    object <- qc.stats(object)
    object <- cc(object)
    object@raw.data <- DATAraw
    object@stats$Phase <- FinalData$assignment.annotation
    return(object)
  }
  ######
  if (return.stats == TRUE) {
    TO <- length(scoring.List)+3
    FinalData <- FinalData[1:TO]
    return(FinalData)
  }
}
