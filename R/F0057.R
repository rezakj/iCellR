#' Demultiplexing HTOs
#'
#' Demultiplexing HTOs
#' @param hto.data HTO raw data
#' @param cov.thr A number which avrage coverage is devided by to set a thershold for low coverage, default = 10.
#' @param assignment.thr A percent above which you decide to set as a good sample assignment/HTO, default = 80.
#' @return An object of class iCellR
#' @examples
#' my.hto <- read.table(file = system.file('extdata', 'dense_umis.tsv',
#'          package = 'iCellR'),
#'          as.is = TRUE)
#' head(my.hto)[1:5]
#'
#' htos <- hto.anno(hto.data = my.hto)
#' head(htos)
#'
#' boxplot(htos$percent.match)
#'
#' @export
hto.anno <- function (hto.data = "data.frame",
                      cov.thr = 10,
                      assignment.thr = 80) {
  # read VDJ data
  data <- hto.data
  data = as.data.frame(t(data))
  mySum <- as.numeric(rowSums(data))
  Neg = mean(mySum)/cov.thr
  myCov = mySum < Neg
  myMax <- as.numeric(apply(data, 1, function(data) {max(data)}))
  myRatio <- (myMax)/(mySum)
  myRatio <- myRatio * 100
  mymat <- myRatio
  mymat[mymat > assignment.thr] <- "good.assignment"
  mymat[mymat < assignment.thr ] <- "unsure"
  myWhich <- as.numeric(apply(data, 1, function(data) {which.max(data)}))
  MyHTOnames = colnames(data)
  myWhich = MyHTOnames[myWhich]
  data <- cbind(data, assignment.annotation = myWhich, percent.match = myRatio, coverage = mySum, low.cov = myCov, assignment.threshold = mymat)
  return(data)
}
