#' Prepare VDJ data
#'
#' This function takes a data frame of VDJ data per cell and prepares it to adds it to the iCellR object.
#' @param vdj.data A data frame containing vdj information.
#' @param cond.name Conditions.
#' @return An object of class iCellR
#' @export
prep.vdj <- function (vdj.data = "data.frame", cond.name = "NULL") {
# read VDJ data
  my.vdj <- vdj.data
  # my.vdj <- subset(my.vdj, productive == "True")
  my.vdj <- subset(my.vdj, raw_clonotype_id != "None")
  cell.barcodes <- my.vdj$barcode
  cell.barcodes <- gsub("-",".",cell.barcodes)
  my.vdj$barcode <- cell.barcodes
  mysum <- (dim(my.vdj)[1]) / 2
  myFreq <- as.data.frame(table(my.vdj$raw_clonotype_id))
  ####
  FreqOrd <- myFreq$Freq
  myFreqMat <- as.matrix(myFreq)
  myFreqMat <- myFreqMat[order(FreqOrd, decreasing = TRUE),]
  myFreq <- as.data.frame(myFreqMat)
  myFreq$Freq <- as.numeric(myFreq$Freq)
  ####
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
