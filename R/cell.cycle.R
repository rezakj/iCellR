#' Cell cycle phase prediction
#'
#' This function takes an object of class iCellR and assignes cell cycle stage for the cells.
#' @param object A data frame containing gene counts for cells.
#' @param scoring.List Genes that are used as a marker for phases.
#' @param return.stats Return the data or object. If FALSE the object would be returned.
#' @param scoring.method Choose from "coverage" or "tirosh" for scoring method. See: https://science.sciencemag.org/content/352/6282/189
#' @return The data frame object
#' @importFrom Hmisc cut2
#' @export
cell.cycle <- function (object = NULL,
                scoring.List = NULL,
                return.stats = FALSE,
                scoring.method = "coverage") {
  if ("iCellR" != class(object)[1]) {
    stop("object should be an object of class iCellR")
  }
  if (is.null(scoring.List)) {
    stop('At least 2 set of genes are required for scoring.List
         Example: scoring.List = c("G0","G1S","G2M","M","MG1","S")')
  }
  ##### get data
  all.genes <- row.names(object@raw.data)
  DATA <- object@raw.data
  object <- cc(object)
#########
if (scoring.method == "coverage") {
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
    dim(subDATA)
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
  if (scoring.method == "tirosh") {
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
    object@stats$Phase <- FinalData$assignment.annotation
    return(object)
  }
  ######
  if (return.stats == TRUE) {
    TO <- length(scoring.List)+1
    FinalData <- FinalData[1:TO]
    return(FinalData)
  }
}
