#' Merge multiple data frames and add the condition names to their cell ids
#'
#' This function takes data frame and merges them while also adding condition names to cell ids..
#' @param samples A character vector of data.frame object names.
#' @param condition.names A character vector of data.frame condition names.
#' @return An object of class iCellR
#' @examples
#' demo <- read.table(
#'         file = system.file('extdata', 'demo_data.txt', package = 'iCellR'),
#'         as.is = TRUE)
#'
#' # Lets divide your sample in to 3 samples as if you have 3 samples and want to merge them.
#' sample1 <- demo[1:30]
#' sample2 <- demo[31:60]
#' sample3 <- demo[61:90]
#'
#' # merge all 3 data and add condition names
#' demo <- data.aggregation(samples =
#'         c("sample1","sample2","sample3"),
#'         condition.names = c("WT","ctrl","KO"))
#' head(demo)[1:4]
#'
#' # make iCellR object
#' myDemo.obj <- make.obj(demo)
#' @importFrom data.table setDT data.table
#' @export
data.aggregation <- function (samples = NULL,
                              condition.names = NULL) {
# check conditions
  if (length(samples) < 2) {
    stop("you need at least 2 samples to merge")
  }
if (length(samples) != length(condition.names)) {
  stop("samples and condition.names need to be of equal size")
}
# fix condition names in the header
  for(i in 1:length(samples)){
    sampledata = get(samples[i])
    colnames(sampledata) <- gsub("_",".",colnames(sampledata))
    colnames(sampledata) <- gsub("-",".",colnames(sampledata))
    colnames(sampledata) <- gsub("/",".or.",colnames(sampledata))
    colnames(sampledata) <- gsub(" ",".",colnames(sampledata))
    colnames(sampledata) <- make.names(colnames(sampledata), unique=TRUE)
    conds = condition.names[i]
    colnames(sampledata) <- paste(conds, colnames(sampledata), sep = "_")
    MyName <- samples[i]
    sampledata <- setDT(sampledata, keep.rownames = TRUE)[]
    eval(call("<-", as.name(MyName), sampledata))
  }
  # merge all by row names
  datalist <- mget(samples)
  mymrgd <- as.data.frame(Reduce(function(x,y) {merge(x,y)}, datalist))
  row.names(mymrgd) <- mymrgd$rn
  # mymrgd <- Reduce(merge, lapply(datalist, function(x) data.frame(x, Row.names = row.names(x))))
  # row.names(mymrgd) <- mymrgd$Row.names
  mymrgd <- mymrgd[,-1]
#  head(mymrgd)[1:5]
      return(mymrgd)
}
