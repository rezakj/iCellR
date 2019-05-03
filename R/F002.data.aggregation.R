#' Merge multiple data frames and add the condition names to their cell ids
#'
#' This function takes data frame and merges them while also adding condition names to cell ids..
#' @param samples A character vector of data.frame object names.
#' @param condition.names A character vector of data.frame condition names.
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' my.data <- data.aggregation(samples = c("sample1","sample2","sample3"),
#'                             condition.names = c("WT","KO","Ctrl"))
#' }
#'
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
    conds = condition.names[i]
    colnames(sampledata) <- paste(conds, colnames(sampledata), sep = "_")
    MyName <- samples[i]
    eval(call("<-", as.name(MyName), sampledata))
  }
  # merge all by row names
  datalist <- mget(samples)
  mymrgd <- Reduce(merge, lapply(datalist, function(x) data.frame(x, Row.names = row.names(x))))
  row.names(mymrgd) <- mymrgd$Row.names
  mymrgd <- mymrgd[,-1]
      return(mymrgd)
}
