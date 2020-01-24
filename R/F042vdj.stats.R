#' VDJ stats
#'
#' This function takes a data frame of VDJ info per cell and dose QC.
#' @param my.vdj A data frame containing VDJ data for cells.
#' @return An object of class iCellR
#' @examples
#' my.vdj <- read.csv(file = system.file('extdata', 'all_contig_annotations.csv',
#'           package = 'iCellR'),
#'           as.is = TRUE)
#' head(my.vdj)
#' dim(my.vdj)
#'
#' My.VDJ <- prep.vdj(vdj.data = my.vdj, cond.name = "NULL")
#' head(My.VDJ)
#' dim(My.VDJ)
#'
#' vdj.stats(My.VDJ)
#'
#' @export
vdj.stats <- function (my.vdj = "data.frame") {
  # read VDJ data
  # chin A
  my.vdj.data <- data.frame(my.vdj$clonotype.Freq,my.vdj$raw_clonotype_id,my.vdj$chain,my.vdj$cdr3)
  my.vdj.data <- subset(my.vdj.data, my.vdj.chain == "TRA")
  ha <- unique((my.vdj.data$my.vdj.clonotype.Freq))
  UpQ <- quantile(ha,0.75)
  LoQ <- quantile(ha,0.25)
  my.vdj.data$Quantile = replace(my.vdj.data$my.vdj.clonotype.Freq, my.vdj.data$my.vdj.clonotype.Freq > UpQ,"Q3")
  my.vdj.data$Quantile = replace(my.vdj.data$Quantile, my.vdj.data$my.vdj.clonotype.Freq < LoQ,"Q1")
  my.vdj.data$Quantile = replace(my.vdj.data$Quantile, my.vdj.data$my.vdj.clonotype.Freq >= LoQ & my.vdj.data$my.vdj.clonotype.Freq <= UpQ,"Mid")
  my.vdj.data$Quantile = replace(my.vdj.data$Quantile, my.vdj.data$my.vdj.clonotype.Freq == 1,"singleton")
  my.vdj.data$Quantile = replace(my.vdj.data$Quantile, my.vdj.data$my.vdj.clonotype.Freq == 2,"doubleton")
  ChainA <- my.vdj.data
  ForPie <- as.data.frame(table(my.vdj.data$Quantile))
  ForPie$Percent <- round((ForPie$Freq/sum(ForPie$Freq))*100,digits = 2)
  ForPie$Name <- paste(ForPie$Var1," (",ForPie$Percent,")", sep = "")
  ForPie$chain <- "TRA"
  # chain B
  my.vdj.data <- data.frame(my.vdj$clonotype.Freq,my.vdj$raw_clonotype_id,my.vdj$chain,my.vdj$cdr3)
  my.vdj.data <- subset(my.vdj.data, my.vdj.chain == "TRB")
  ha <- unique((my.vdj.data$my.vdj.clonotype.Freq))
  UpQ <- quantile(ha,0.75)
  LoQ <- quantile(ha,0.25)
  my.vdj.data$Quantile = replace(my.vdj.data$my.vdj.clonotype.Freq, my.vdj.data$my.vdj.clonotype.Freq > UpQ,"Q3")
  my.vdj.data$Quantile = replace(my.vdj.data$Quantile, my.vdj.data$my.vdj.clonotype.Freq < LoQ,"Q1")
  my.vdj.data$Quantile = replace(my.vdj.data$Quantile, my.vdj.data$my.vdj.clonotype.Freq >= LoQ & my.vdj.data$my.vdj.clonotype.Freq <= UpQ,"Mid")
  my.vdj.data$Quantile = replace(my.vdj.data$Quantile, my.vdj.data$my.vdj.clonotype.Freq == 1,"singleton")
  my.vdj.data$Quantile = replace(my.vdj.data$Quantile, my.vdj.data$my.vdj.clonotype.Freq == 2,"doubleton")
  ChainB <- my.vdj.data
  ForPie1 <- as.data.frame(table(my.vdj.data$Quantile))
  ForPie1$Percent <- round((ForPie1$Freq/sum(ForPie1$Freq))*100,digits = 2)
  ForPie1$Name <- paste(ForPie1$Var1," (",ForPie1$Percent,")", sep = "")
  ForPie1$chain <- "TRB"
  # both
  FullPie <- rbind(ForPie,ForPie1)
  FullPie$Var1 <- factor(FullPie$Var1, levels = c("singleton","doubleton","Q1","Mid","Q3"))
  myPIE <- ggplot(FullPie,aes(y=Freq, x="", fill = Var1, group= chain)) +
    geom_bar(stat = "identity", position = "fill") + theme_bw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + coord_polar(theta="y") + facet_wrap(~ chain, ncol =1)
  #
  data <- (unique(subset(ChainB, ChainB$Quantile == "Q3")))
  data <- (data[order(data$my.vdj.clonotype.Freq, decreasing = TRUE),])
  data1 <- (unique(subset(ChainA, ChainA$Quantile == "Q3")))
  data1 <- (data1[order(data1$my.vdj.clonotype.Freq, decreasing = TRUE),])
  DATA <- rbind(data,data1)
  DATA <- (DATA[order(DATA$my.vdj.clonotype.Freq, decreasing = TRUE),])
  colnames(DATA) <- c("freq","colonotype","chain","cdr3","Q")
  DATA$cdr3 <- factor(DATA$cdr3, levels = rev(unique(DATA$cdr3)))
  freq = log2(DATA$freq)
  myPLOT <- ggplot(DATA,aes(x= freq, y=cdr3, col=chain)) +
    geom_line() + geom_point(size = freq) +
    theme_bw(base_size = 16)
  #  write.table((my.vdj),file="my.vdj.tsv",sep="\t",row.names =F)
  #  grid.arrange(myPLOT,myPIE, ncol=2, widths=c(2,1), heights=c(10,1))
  return(grid.arrange(myPLOT,myPIE, ncol=2))
}
