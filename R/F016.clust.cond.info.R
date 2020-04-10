#' Calculate cluster and conditions frequencies
#'
#' This function takes an object of class iCellR and calculates cluster and conditions frequencies.
#' @param x An object of class iCellR.
#' @param plot.type Choose from pie/pie.cond or bar/bar.cond, defult = pie.
#' @param my.out.put Chose from "data" or "plot", default = "data".
#' @param normalize.ncell If TRUE the values will be normalized to the number of cells by downsampling.
#' @param normalize.by Chose from "sf" (size factor) or "percentage", default = "percentage".
#' @return An object of class iCellR.
#' @export
clust.cond.info <- function (x = NULL,
                             plot.type = "pie",
                             my.out.put = "data",
                             normalize.ncell = TRUE,
                             normalize.by = "percentage") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ###################
  Cells <- colnames(x@main.data)
  MYConds <- as.character((unique(data.frame(do.call('rbind', strsplit(as.character(Cells),'_',fixed=TRUE)))[1]))$X1)
  if (length(MYConds) == 1) {
    stop("You need more then one condition/sample to run this function")
  }
  if (length(MYConds) == 0) {
    stop("You need more then one condition/sample to run this function")
  }
  ###################
  DATA <- (x@best.clust)
  Conds <- (as.data.frame(do.call("rbind", strsplit(row.names(DATA), "_")))[1])
  ForNorm1 <- as.data.frame(table(Conds))
  ForNorm <- min(ForNorm1$Freq)
  SizeFactors <- round(ForNorm1$Freq/ForNorm,3)
  ForNorm1$SF <- SizeFactors
  #
  clusts <- (as.data.frame(DATA$clusters))
  cond.clust <- cbind(Conds, clusts)
  colnames(cond.clust) <- c("conditions","clusters")
  Conds <- as.character(ForNorm1$Conds)
  My.Conds.data <- cond.clust
  #################
  DATA <- as.data.frame(table(cond.clust))
  Freq <- DATA$Freq
  colnames(ForNorm1) <- c("conditions","TC","SF")
  DATA <- merge(ForNorm1,DATA,by="conditions")
  DATA$Norm.Freq <- round(DATA$Freq/DATA$SF,3)
  DATA$percentage <- round((DATA$Freq/DATA$TC)*100,2)
################
# as.data.frame(table(Conds))
  # bar
  myBP <- ggplot(DATA,aes(y=Freq, x=conditions, fill = conditions)) +
    geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x=element_text(angle=90)) + facet_wrap(~ clusters, scales = "free")
  # bar2
  myBP2 <- ggplot(DATA,aes(y=Freq, x=conditions, fill = clusters)) +
    geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x=element_text(angle=90)) + facet_wrap(~ conditions, scales = "free")
    # pie
  myPIE <- ggplot(DATA,aes(y=Freq, x="", fill = conditions)) +
    geom_bar(stat = "identity", position = "fill") + theme_bw() + facet_wrap(~ clusters) +
      theme(axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) + coord_polar(theta="y")
  # pie2
  myPIE2 <- ggplot(DATA,aes(y=Freq, x="", fill = clusters)) +
    geom_bar(stat = "identity", position = "fill") + theme_bw() + facet_wrap(~ conditions) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + coord_polar(theta="y")

  #############
  if (normalize.ncell == TRUE) {
    if(normalize.by == "sf") {
    # bar
    myBP <- ggplot(DATA,aes(y=Norm.Freq, x=conditions, fill = conditions)) +
      geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x=element_text(angle=90)) + facet_wrap(~ clusters, scales = "free")
    # bar2
    myBP2 <- ggplot(DATA,aes(y=Norm.Freq, x=conditions, fill = clusters)) +
      geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x=element_text(angle=90)) + facet_wrap(~ conditions, scales = "free")
    # pie
    myPIE <- ggplot(DATA,aes(y=Norm.Freq, x="", fill = conditions)) +
      geom_bar(stat = "identity", position = "fill") + theme_bw() + facet_wrap(~ clusters) +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) + coord_polar(theta="y")
    # pie2
    myPIE2 <- ggplot(DATA,aes(y=Norm.Freq, x="", fill = clusters)) +
      geom_bar(stat = "identity", position = "fill") + theme_bw() + facet_wrap(~ conditions) +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) + coord_polar(theta="y")
    }
    if (normalize.by == "percentage") {
      # bar
      myBP <- ggplot(DATA,aes(y=percentage, x=conditions, fill = conditions)) +
        geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x=element_text(angle=90)) + facet_wrap(~ clusters, scales = "free")
      # bar2
      myBP2 <- ggplot(DATA,aes(y=percentage, x=conditions, fill = clusters)) +
        geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x=element_text(angle=90)) + facet_wrap(~ conditions, scales = "free")
      # pie
      myPIE <- ggplot(DATA,aes(y=percentage, x="", fill = conditions)) +
        geom_bar(stat = "identity", position = "fill") + theme_bw() + facet_wrap(~ clusters) +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) + coord_polar(theta="y")
      # pie2
      myPIE2 <- ggplot(DATA,aes(y=percentage, x="", fill = clusters)) +
        geom_bar(stat = "identity", position = "fill") + theme_bw() + facet_wrap(~ conditions) +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) + coord_polar(theta="y")
    }
  }
#  write.table((DATA), file="clust_cond_freq_info.txt", sep="\t", row.names =FALSE)
#  message("clust_cond_freq_info.txt file has beed generated.")
  if (my.out.put == "plot") {
    if (plot.type == "bar") {
      return(myBP)
    }
    if (plot.type == "bar.cond") {
      return(myBP2)
    }
    if (plot.type == "pie") {
      return(myPIE)
    }
    if (plot.type == "pie.cond") {
      return(myPIE2)
    }
  }
#
  if (my.out.put == "data") {
    attributes(x)$my.freq <- DATA
    return(x)
  }
}
