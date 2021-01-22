#' Add image data to iCellR object
#'
#' This function takes a list of image data adds it to the iCellR object.
#' @param x An object of class iCellR.
#' @param image.data.list A character vector of list object names. Lists should be made using "image.capture.10x" function. .
#' @param condition.names A character vector of condition names.
#' @return An object of class iCellR
#' @export
add.10x.image <- function (x = NULL,
                           image.data.list = NULL,
                           condition.names = NULL) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  #####
  if (length(image.data.list) == 1) {
    message(" ... A singe image sample is added to the iCellR object.")
    if (!is.null(condition.names)) {
      stop("condition names can't be given for a single sample")
    }
    data <- get(image.data.list)
    attributes(x)$spatial.info <- data
    data <- data$position
    rownames(data) <- gsub("-",".",rownames(data))
    rownames(data) <- gsub("_",".",rownames(data))
    attributes(x)$spatial.data <- data
    return(x)
    }
########## more than one conditions
  if (length(image.data.list) > 1) {
    ToPrint2 <- paste("    ...", length(image.data.list),"samples are added to the iCellR object.", sep=" ")
    message(ToPrint2)
    if (!is.null(condition.names)) {
      ToPrint2 <- paste(as.character(condition.names), collapse=",")
      ToPrint2 <- paste("         condition names are:", ToPrint2, sep=" ")
      message(ToPrint2)
    }
    if (length(image.data.list) != length(condition.names)) {
      stop("condition names and number of lists don't match")
    }
    for(i in 1:length(image.data.list)){
      data = get(image.data.list[i])
      data <- data$position
      rownames(data) <- gsub("-",".",rownames(data))
      rownames(data) <- gsub("_",".",rownames(data))
      conds = condition.names[i]
      rownames(data) <- paste(conds, rownames(data), sep = "_")
      MyName <- paste("SampleName",conds, sep="_")
      eval(call("<-", as.name(MyName), data))
    }
##########
    filenames <- ls(pattern="SampleName_")
    datalist <- mget(filenames)
    names(datalist) <- NULL
    data <- do.call("rbind", datalist)
    attributes(x)$spatial.data <- data
############
    for(i in 1:length(image.data.list)){
      data = get(image.data.list[i])
      conds = condition.names[i]
      MyName <- paste("SampleName",conds, sep="_")
      eval(call("<-", as.name(MyName), data))
    }
    filenames <- ls(pattern="SampleName_")
    datalist <- mget(filenames)
    attributes(x)$spatial.info <- datalist
    return(x)
  }
}
