#' Add CITE-seq antibody-derived tags (ADT)
#'
#' This function takes a data frame of ADT values per cell and adds it to the iCellR object.
#' @param x An object of class iCellR.
#' @param adt.data A data frame containing ADT counts for cells.
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' my.obj <- add.adt(my.obj, adt.data = adt.data)
#' }
#'
#' @export
prep.vdj <- function (vdj.data = "all_contig_annotations.csv", cond.name = "NULL") {
# read VDJ data
  my.vdj <- read.table(vdj.data, header = T, sep=",")
  my.vdj <- subset(my.vdj, productive == "True")
  my.vdj <- subset(my.vdj, raw_clonotype_id != "None")
  mysum <- (dim(my.vdj)[1]) / 2
  myFreq <- as.data.frame(table(my.vdj$raw_clonotype_id))
  myFreq <- myFreq[order(myFreq$Freq, decreasing = T),]
  myFreq$Freq <- (round(myFreq$Freq / 2))
  colono.sum <- dim(myFreq)[1]
  myFreq$my.raw_clonotype_id <- myFreq$Var1
  myFreq$clonotype.Freq <- myFreq$Freq
  myFreq$proportion <- (myFreq$Freq / mysum)
  myFreq$total.colonotype <- colono.sum
  myFreq$Freq <- NULL
  my.vdj <- merge(my.vdj,myFreq, by.x="raw_clonotype_id",by.y="Var1", all.x=T, all.y=F)
  if (cond.name != "NULL") {
    my.vdj$barcode <- paste(cond.name, my.vdj$barcode, sep = "_")
  }
#  write.table((my.vdj),file="my.vdj.tsv",sep="\t",row.names =F)
  return(my.vdj)
}
