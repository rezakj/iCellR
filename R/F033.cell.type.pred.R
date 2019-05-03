#' Create heatmaps or dot plots for genes in clusters to find thier cell types using ImmGen data.
#'
#' This function takes an object of class iCellR and genes and provides a heatmap.
#' @param immgen.data Choose from "rna", "uli.rna" or "mca", default = "rna"
#' @param gene A set of gene names to used to predict cell type.
#' @param plot.type Choose from "heatmap" od "point.plot", default = "heatmap"
#' @param top.cell.types Top cell types sorted by cumulative expression, default = 25.
#' @param heat.colors Colors for heatmap, default = c("blue" ,"white", "red").
#' @return An object of class iCellR
#' @examples
#' \dontrun{
#' imm.gen(immgen.data = "uli.rna", gene = MyGenes, plot.type = "heatmap")
#' imm.gen(immgen.data = "rna", gene = MyGenes, plot.type = "point.plot")
#' }
#' @import pheatmap
#' @export
cell.type.pred <- function (immgen.data = "rna",
                     gene = "NULL",
                     top.cell.types = 50,
                     plot.type = "heatmap",
                     heat.colors = c("blue","white", "red")) {
  ## get main data
if (immgen.data == "rna") {
  my.data <- immgen.rna
  MyTitle = "Top 23 out of 23 (ImmGen RNA-seq)"
  }
  if (immgen.data == "uli.rna") {
    my.data <- immgen.uli.rna
    MyTitle = paste("Top", top.cell.types, "out of 157 (ImmGen ultra-low-input RNA-seq)")
  }
  if (immgen.data == "mca") {
    my.data <- mouse.cell.atlas
    MyTitle = paste("Top", top.cell.types, "out of 157 (ImmGen ultra-low-input RNA-seq)")
  }
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
  MYdf <- as.data.frame(sort(colSums(my.data),decreasing = F))
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
                    col = mycol,
                    show_colnames = T,
                    cluster_rows = T,
                    cluster_cols = T,
                    scale = "row"))
  }
}
