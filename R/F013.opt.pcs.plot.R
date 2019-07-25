#' Find optimal number of PCs for clustering
#'
#' This function takes an object of class iCellR and finds optimal number of PCs for clustering.
#' @param x An object of class iCellR.
#' @param pcs.in.plot Number of PCs to show in plot, defult = 50.
#' @return An object of class iCellR.
#' @examples
#' opt.pcs.plot(demo.obj)
#'
#' @export
opt.pcs.plot <- function (x = NULL, pcs.in.plot = 50) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }

#  DATA <- x@pca.data
  DATA <- x@pca.info$sdev
  OPTpcs <- mean(DATA)*2
  OPTpcs <- (DATA > OPTpcs)
  OPTpcs <- length(OPTpcs[OPTpcs==TRUE]) + 1
  Titel <- paste("Optimal number of PCs (1:", OPTpcs, ")", sep="")
  # fix DATA
  DATA <- as.data.frame(DATA)
  w = log1p(1:100)
  Rows <- c(1:dim(DATA)[1])
  DATA <- cbind(Rows, DATA)
  colnames(DATA) <- c("PCs","SDs")
  DATA <- head(DATA,pcs.in.plot)
### plot
  LogSDs <- log2(DATA$SDs)
  Log2PCs <- log2(DATA$PCs)
  myPLOT <- ggplot(DATA,aes(y=SDs,x=PCs, col = LogSDs)) +
    geom_line() + geom_point(size = LogSDs) +
    scale_x_continuous(trans = "log1p") +
    scale_y_continuous(trans = "log1p") +
    geom_vline(xintercept = OPTpcs,linetype="dotted") +
    theme_bw(base_size = 16)+
    xlab("Principal Components") +
    ylab("PC Standard Deviations") +
    scale_colour_gradient(low = "gray", high = "red", name="log2 (SDs)") +
    ggtitle(Titel)
# return
  return(myPLOT)
}

