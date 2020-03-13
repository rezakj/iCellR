#' Create heatmaps or dot plots for genes in clusters to find thier cell types using ImmGen data.
#'
#' This function takes an object of class iCellR and genes and provides a heatmap.
#' @param immgen.data Choose from ,"GSE109125","GSE122108","GSE122597","GSE124829","GSE15907","GSE37448", rna", "uli.rna" or "mca", default = "rna"
#' @param gene A set of gene names to used to predict cell type.
#' @param plot.type Choose from "heatmap" od "point.plot", default = "heatmap"
#' @param top.cell.types Top cell types sorted by cumulative expression, default = 25.
#' @param heat.colors Colors for heatmap, default = c("blue" ,"white", "red").
#' @return An object of class iCellR
#' @import pheatmap
#' @export
cell.type.pred <- function (immgen.data = "rna",
                     gene = "NULL",
                     top.cell.types = 50,
                     plot.type = "heatmap",
                     heat.colors = c("blue","white", "red")) {
  ## get main data
  if (immgen.data == "GSE109125") {
    if (!exists("Immgen.GSE109125.205")) {
      stop("Please download the 'Immgen.GSE109125.205.rda' file from 'https://genome.med.nyu.edu/results/external/iCellR/data/' and load it like: load('Immgen.GSE109125.205.rda')")
    }
    my.data <- Immgen.GSE109125.205
    MyTitle = paste("Top", top.cell.types, "out of 205 (ImmGen RNA-seq)")
  }
  ###
  if (immgen.data == "GSE122108") {
    if (!exists("Immgen.GSE122108.412")) {
      stop("Please download the 'Immgen.GSE122108.412.rda' file from 'https://genome.med.nyu.edu/results/external/iCellR/data/' and load it like: load('Immgen.GSE122108.412.rda')")
    }
    my.data <- Immgen.GSE122108.412
    MyTitle = paste("Top", top.cell.types, "out of 412 (ImmGen RNA-seq)")
  }
  ###
  if (immgen.data == "GSE122597") {
    if (!exists("Immgen.GSE122597.83")) {
      stop("Please download the 'Immgen.GSE122597.83.rda' file from 'https://genome.med.nyu.edu/results/external/iCellR/data/' and load it like: load('Immgen.GSE122597.83.rda')")
    }
    my.data <- Immgen.GSE122597.83
    MyTitle = paste("Top", top.cell.types, "out of 83 (ImmGen RNA-seq)")
  }
  ###
  if (immgen.data == "GSE124829") {
    if (!exists("Immgen.GSE124829.190")) {
      stop("Please download the 'Immgen.GSE124829.190.rda' file from 'https://genome.med.nyu.edu/results/external/iCellR/data/' and load it like: load('Immgen.GSE124829.190.rda')")
    }
    my.data <- Immgen.GSE124829.190
    MyTitle = paste("Top", top.cell.types, "out of 190 (ImmGen RNA-seq)")
  }
  ###
  if (immgen.data == "GSE37448") {
    if (!exists("Immgen.microarray.GSE37448.189")) {
      stop("Please download the 'Immgen.microarray.GSE37448.189.rda' file from 'https://genome.med.nyu.edu/results/external/iCellR/data/' and load it like: load('Immgen.microarray.GSE37448.189.rda')")
    }
    my.data <- Immgen.microarray.GSE37448.189
    MyTitle = paste("Top", top.cell.types, "out of 189 (ImmGen Microarray RNA-seq)")
  }
  ###
  if (immgen.data == "GSE15907") {
    if (!exists("Immgen.microarray.GSE15907.653")) {
      stop("Please download the 'Immgen.microarray.GSE15907.653.rda' file from 'https://genome.med.nyu.edu/results/external/iCellR/data/' and load it like: load('Immgen.microarray.GSE15907.653.rda')")
    }
    my.data <- Immgen.microarray.GSE15907.653
    MyTitle = paste("Top", top.cell.types, "out of 653 (ImmGen Microarray RNA-seq)")
  }
  ####
if (immgen.data == "rna") {
  if (!exists("immgen.rna")) {
    stop("Please download the 'immgen.rna.rda' file from 'https://genome.med.nyu.edu/results/external/iCellR/data/' and load it like: load('immgen.rna.rda')")
  }
  my.data <- immgen.rna
  MyTitle = "Top 23 out of 23 (ImmGen RNA-seq)"
  }
  if (immgen.data == "uli.rna") {
    if (!exists("immgen.uli.rna")) {
      stop("Please download the 'immgen.uli.rna.rda' file from 'https://genome.med.nyu.edu/results/external/iCellR/data/' and load it like: load('immgen.uli.rna.rda')")
    }
    my.data <- immgen.uli.rna
    MyTitle = paste("Top", top.cell.types, "out of 157 (ImmGen ultra-low-input RNA-seq)")
  }
  if (immgen.data == "mca") {
    if (!exists("mouse.cell.atlas")) {
      stop("Please download the 'mouse.cell.atlas.rda' file from 'https://genome.med.nyu.edu/results/external/iCellR/data/' and load it like: load('mouse.cell.atlas.rda')")
    }
    my.data <- mouse.cell.atlas
    MyTitle = paste("Top", top.cell.types, "Mouse Cell Atlas")
  }
  ########################################
  # fix row names
  MyRows <- toupper(row.names(my.data))
  MyRows <- gsub("-",".", MyRows)
  rownames(my.data) <- MyRows
# get genes
  gene <- toupper(gene)
  my.data <- subset(my.data,row.names(my.data) %in% gene)
  # scale
  my.data <- log2(my.data + 1)
  # gg plot
  MYdf <- as.data.frame(sort(colSums(my.data),decreasing = FALSE))
  colnames(MYdf) <- c("MyLev")
  MYdf <- cbind(cells = rownames(MYdf),MYdf)
  MYdf <- tail(MYdf, top.cell.types)
  MYdf$cells <- factor(MYdf$cells, levels = row.names(MYdf))
  #
  if (plot.type == "point.plot") {
#    mySize = log2(MYdf$MyLev)
    return(ggplot(MYdf, aes(x = log2(MyLev),y=cells, col = log2(MyLev), size =log2(MyLev))) +
             geom_point() +
             xlab("log2(cumulative expression)") +
             ylab("ImmGen RNA-seq Cell Types") +
             ggtitle(MyTitle) +
             theme_classic() +
             scale_size("") +
             scale_colour_gradient(low = "gray", high = "red", name="")
             )
  }
  # heatmap colors
  mycol <- colorRampPalette(heat.colors)(n = 100)
  # return
  if (plot.type == "heatmap") {
    return(pheatmap(as.matrix(my.data),
                    color = mycol,
                    show_colnames = TRUE,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    scale = "row"))
  }
}
