#' Run PCA on the main data
#'
#' This function takes an object of class iCellR and runs PCA on the main data.
#' @param x An object of class iCellR.
#' @param clust.method Choose from "base.mean.rank" or "gene.model", defult is "base.mean.rank".
#' @param top.rank A number taking the top genes ranked by base mean, defult = 500.
#' @param gene.list A list of genes to be used for PCA. If "clust.method" is set to "gene.model", defult = "my_model_genes.txt".
#' @return An object of class iCellR.
#' @examples
#' \dontrun{
#' my.obj <- run.pca(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")
#' }
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
    data <- data[order(data[1], decreasing = T),]
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
  write.table((best.genes),file="my_model_PC_genes.txt", row.names =F, quote = F, col.names = F)
  print("my_model_PC_genes.txt file is generated, which can be used for clustering.")
}
