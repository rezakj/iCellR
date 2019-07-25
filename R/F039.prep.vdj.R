#' Prepare VDJ data
#'
#' This function takes a data frame of VDJ data per cell and prepares it to adds it to the iCellR object.
#' @param vdj.data A data frame containing vdj information.
#' @param cond.name Conditions.
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' Read your VDJ data (in this case in VDJ.tsv file) and add to your object as below
#'
#' my.vdj.data <- read.table("VDJ.tsv")
#'
#' VDJ <- prep.vdj(my.obj, adt.data = my.vdj.data)
#'
#' head(VDJ)
#' }
#'
#' @export
prep.vdj <- function (vdj.data = "all_contig_annotations.csv", cond.name = "NULL") {
# read VDJ data
  my.vdj <- read.table(vdj.data, header = TRUE, sep=",")
  my.vdj <- subset(my.vdj, productive == "True")
  my.vdj <- subset(my.vdj, raw_clonotype_id != "None")
  mysum <- (dim(my.vdj)[1]) / 2
  myFreq <- as.data.frame(table(my.vdj$raw_clonotype_id))
  myFreq <- myFreq[order(myFreq$Freq, decreasing = TRUE),]
  myFreq$Freq <- (round(myFreq$Freq / 2))
  colono.sum <- dim(myFreq)[1]
  myFreq$my.raw_clonotype_id <- myFreq$Var1
  myFreq$clonotype.Freq <- myFreq$Freq
  myFreq$proportion <- (myFreq$Freq / mysum)
  myFreq$total.colonotype <- colono.sum
  myFreq$Freq <- NULL
  my.vdj <- merge(my.vdj,myFreq, by.x="raw_clonotype_id",by.y="Var1", all.x=TRUE, all.y=FALSE)
  if (cond.name != "NULL") {
    my.vdj$barcode <- paste(cond.name, my.vdj$barcode, sep = "_")
  }
#  write.table((my.vdj),file="my.vdj.tsv",sep="\t",row.names =F)
  return(my.vdj)
}
