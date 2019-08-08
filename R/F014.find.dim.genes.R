#' Find model genes from PCA data
#'
#' This function takes an object of class iCellR finds the model genes to run a second round of PCA.
#' @param x An object of class iCellR.
#' @param dims PC dimentions to be used.
#' @param top.pos Number of top positive marker genes to be taken from each PC, default = 15.
#' @param top.neg Number of top negative marker genes to be taken from each PC, default = 5.
#' @return An object of class iCellR.
#' @examples
#'
#' demo.obj <- find.dim.genes(demo.obj, dims = 1:10,top.pos = 20, top.neg = 20)
#'
#' head(demo.obj@gene.model)
#'
#' @export
find.dim.genes <- function (x = NULL,
                            dims = 1:10,
                            top.pos = 15,
                            top.neg = 5) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  # geth the genes and scale them based on model
  DATA <- as.data.frame(x@pca.info$x)
  DATA <- DATA[dims]
  PCs <- dim(DATA)[2]
  PCs <- c(1:PCs)
  ##
  for (i in PCs) {
    data <- DATA[i]
    data <- cbind(data, genes=row.names(data))
    data <- data[order(data[1], decreasing = TRUE),]
    TOP <- row.names(head(data,top.pos))
    BOT <- row.names(tail(data,top.neg))
    MYgenes <- c(TOP,BOT)
    NameCol=paste("PCgenes",i,sep="_")
    eval(call("<-", as.name(NameCol), MYgenes))
  }
#
  filenames <- ls(pattern="PCgenes_")
  datalist <- mget(filenames)
  best.genes <- unique(as.character(as.matrix(as.data.frame(datalist))))
#
#  write.table((best.genes),file="my_model_PC_genes.txt", row.names =FALSE, quote = FALSE, col.names = FALSE)
#  message("my_model_PC_genes.txt file is generated, which can be used for clustering.")
  attributes(x)$gene.model <- as.character(as.matrix(best.genes))
  return(x)
  }
