#' iCellR Batch Alignment (IBA)
#'
#' This function takes an object of class iCellR and runs CCCA or CPCA batch alignment.
#' @param x An object of class iCellR.
#' @param ba.method Batch alignment method. Choose from "CCCA" and "CPCA", default = "CPCA".
#' @param dims PC dimentions to be used
#' @param k number of neighboring cells for KNN, default = 10.
#' @param method Choose from "base.mean.rank" or "gene.model", default is "base.mean.rank". If gene.model is chosen you need to provide gene.list.
#' @param top.rank A number. Taking the top genes ranked by base mean, default = 500.
#' @param plus.log.value A number to add to each value in the matrix before log transformasion to aviond Inf numbers, default = 0.1.
#' @param gene.list A charactor vector of genes to be used for PCA. If "clust.method" is set to "gene.model", default = "my_model_genes.txt".
#' @param scale.data If TRUE the data will be scaled (log2 + plus.log.value), default = TRUE.
#' @import progress
#' @export
iba <- function (x = NULL,
                 dims = 1:30,
                 k = 10,
                 ba.method = "CPCA",
                 method = "base.mean.rank",
                 top.rank = 500,
                 plus.log.value = 0.1,
                 scale.data = TRUE,
                 gene.list = "character") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  if (ba.method != "CCCA" & ba.method != "CPCA") {
    stop("method should be CPCA or CCCA")
  }
  # get data
  start_time1 <- Sys.time()
  DATA <- x@main.data
  my.data = DATA
  myNN = k
  message(paste(" main data dimentions:",dim(DATA)[1],"genes and",dim(DATA)[2],"cells"))
  #########
  #########################
  start_time <- Sys.time()
  message(paste("   Running PCA ..."))
  x <- run.pca(x,
               method = method,
               top.rank = top.rank,
               plus.log.value = plus.log.value,
               scale.data = scale.data,
               gene.list = gene.list)
  end_time <- Sys.time()
  Time = difftime(end_time,start_time,units = "mins")
  Time = round(as.numeric(Time),digits = 2)
  message(paste("   PCA performed in",Time,"mins"))
  num.of.PCs = c(dims)
  my.data.my.pca = t(x@pca.info$rotation)[num.of.PCs, ]
  #########
  start_time <- Sys.time()
  message(paste("   Calculating distance ..."))
  My.distances = as.matrix(dist(t(my.data.my.pca)))
  end_time <- Sys.time()
  Time = difftime(end_time,start_time,units = "mins")
  Time = round(as.numeric(Time),digits = 2)
  message(paste("   Calculated distance in",Time,"mins"))
  ncells = dim(my.data)[2]
  ########################
  ### get conditions
  ha <- colnames(My.distances)
  ha <- data.frame(do.call('rbind', strsplit(as.character(ha),'_',fixed=TRUE)))[1]
  ha <- unique(as.character(as.matrix(ha)))
  conditions <- ha
  #####
  #### function
  message(paste("   Calculating joint distance ... "))
  ### time
  pb <- progress_bar$new(total = ncells,
                         format = "[:bar] :current/:total (:percent) :elapsedfull eta: :eta",
                         clear = FALSE, width= 60)
  ### end time
  dataha <- as.data.frame(row.names(My.distances))
  colnames(dataha) <- "Ids"
  ####
  KNN1 = lapply(1:ncells, function(findKNN){
    ###### time
    pb$tick()
    #### end time
    ######### loop
    for(i in conditions){
      ha <- paste("^",i,"_",sep="")
#      CellOrd <- colnames(My.distances)[(GETord(My.distances[,findKNN]))]
      CellOrd <- colnames(My.distances)[(order(My.distances[,findKNN]))]
      CellsId <- subset(CellOrd, grepl(ha, CellOrd))[1:myNN]
      #          CellsId <- grep(ha,CellOrd,value=T, invert=F)[1:myNN]
      NameCol=paste("MySet",i,sep="_")
      eval(call("<-", as.name(NameCol), CellsId))
    }
    filenames <- ls(pattern="MySet_")
    do <- as.character(unlist(mget(filenames)))
    ########## end loop
    #do <- paste("^",do,"$",sep="")
    #do <- paste(do,collapse="|")
    #do <- grep(do,row.names(My.distances),value=F)
    do <- as.numeric(rownames(subset(dataha,dataha$Ids %in% do)))
  })
  #### end of function
  message(paste("  "))
  ######
  if (ba.method == "CPCA") {
    message(paste("   Alignment method: Combined Principal Component Alignment (CPCA)"))
    my.data <- as.data.frame(my.data.my.pca)
  }
  if (ba.method == "CCCA") {
    message(paste("   Alignment method: Combined Coverage Correction Alignment (CCCA)"))
  }
  ######
  ### time
  my.data <- as.matrix(my.data)
  pb <- progress_bar$new(total = ncells,
                         format = "[:bar] :current/:total (:percent) :elapsedfull eta: :eta",
                         clear = FALSE, width= 60)
  ### end time
  data.sum1 = sapply(KNN1, function(sum.cov){
    pb$tick()
#    GETmean(my.data[, sum.cov])})
    rowMeans(my.data[, sum.cov])})
  ############
  data.sum1 <- as.data.frame(data.sum1)
  row.names(data.sum1) <- row.names(my.data)
  colnames(data.sum1) <- colnames(my.data)
#  data.sum1 <- round(data.sum1, digits = 3)
  #### order samples (order colnames by the original data)
  BestOrder <- colnames(DATA)
  data.sum1 <- data.sum1[ , order(match(names(data.sum1),BestOrder))]
  ########
  attributes(x)$imputed.data <- data.sum1
  start_time <- Sys.time()
  message(paste(" PCA ..."))
  if (ba.method == "CPCA") {
    x <- run.pca(x,data.type = "imputed", scale.data = FALSE)
  }
  if (ba.method == "CCCA") {
    x <- run.pca(x,data.type = "imputed",
                 method = method,
                 top.rank = top.rank,
                 plus.log.value = plus.log.value,
                 scale.data = scale.data,
                 gene.list = gene.list)
  }
  end_time <- Sys.time()
  Time = difftime(end_time,start_time,units = "mins")
  Time = round(as.numeric(Time),digits = 2)
  message(paste("   PCA performed in",Time,"mins"))
  message(paste("All done!"))
  end_time1 <- Sys.time()
  Time = difftime(end_time1,start_time1,units = "mins")
  Time = round(as.numeric(Time),digits = 2)
  message(paste("Total time",Time,"mins"))
  #######
  if (ba.method == "CPCA") {
    attributes(x)$imputed.data <- as.data.frame(NULL)
  }
  return(x)
}
