#' Make a gene model for clustering
#'
#' This function takes an object of class iCellR and provides a gene list for clustering based on the parameters set in the model.
#' @param x An object of class iCellR.
#' @param dispersion.limit A number for taking the genes that have dispersion above this number, default = 1.5.
#' @param base.mean.rank A number taking the top genes ranked by base mean, default = 500.
#' @param gene.num.max Maximum number of genes , default = 2000.
#' @param non.sig.col Color for the genes not used for the model, default = "darkgray".
#' @param right.sig.col Color for the genes above the dispersion limit, default = "chartreuse3".
#' @param left.sig.col Color for the genes above the rank limit, default = "cadetblue3".
#' @param disp.line.col Color of the line for dispersion limit, default = "black".
#' @param rank.line.col Color of the line for rank limit, default = "red".
#' @param cell.size A number for the size of the points in the plot, default = 1.75.
#' @param no.mito.model If set to TRUE, mitochondrial genes would be excluded from the gene list made for clustering, default = TRUE.
#' @param mark.mito Mark mitochondrial genes in the plot, default = TRUE.
#' @param my.out.put Chose from "data" or "plot", default = "data".
#' @param no.cell.cycle If TRUE the cell cycle genes will be removed (s.phase and g2m.phase), default = TRUE.
#' @param cell.transparency Color transparency for the points in the plot, default = 0.5.
#' @param interactive If set to TRUE an interactive HTML file will be created, default = TRUE.
#' @param out.name If "interactive" is set to TRUE, the out put name for HTML, default = "plot".
#' @return An object of class iCellR.
#' @examples
#' make.gene.model(demo.obj,
#'                dispersion.limit = 1.5,
#'                base.mean.rank = 500,
#'                no.mito.model = TRUE,
#'                mark.mito = TRUE,
#'                interactive = FALSE,
#'                my.out.put = "plot",
#'                out.name = "gene.model")
#'
#' demo.obj <- make.gene.model(demo.obj,
#'                            dispersion.limit = 1.5,
#'                            base.mean.rank = 500,
#'                            no.mito.model = TRUE,
#'                            mark.mito = TRUE,
#'                            interactive = FALSE,
#'                            out.name = "gene.model")
#'
#' head(demo.obj@gene.model)
#'
#' @import ggrepel
#' @export
make.gene.model <- function (x = NULL,
                             dispersion.limit = 1.5,
                             base.mean.rank = 500,
                             gene.num.max = 2000,
                             non.sig.col = "darkgray",
                             right.sig.col = "chartreuse3",
                             left.sig.col = "cadetblue3",
                             disp.line.col = "black",
                             rank.line.col = "red",
                             my.out.put = "data",
                             cell.size = 1.75,
                             cell.transparency = 0.5,
                             no.mito.model = TRUE,
                             no.cell.cycle = TRUE,
                             mark.mito = TRUE,
                             interactive = TRUE,
                             out.name = "plot") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  DATA <- x@gene.data
#  x = c(1:100)
#  y = log1p(1:100)
  # get conditions
  Myconds <- as.character(unique(DATA$condition))
  data <- subset(DATA, DATA$condition == "all")
  # get mito genes
  mito.genes <- grep(pattern = "^mt\\-", x = data$genes, value = TRUE, ignore.case = TRUE)
  if ( length(mito.genes) == 0 ) {
    mito.genes <- grep(pattern = "^mt\\.", x = data$genes, value = TRUE, ignore.case = TRUE)
  }
  top_labelled <- subset(data, data$gene %in% mito.genes)
  # get cell cycle genes
  s.phase.genes <- s.phase
  g2m.phase.genes <- g2m.phase
  sg2m <- c(s.phase.genes,g2m.phase.genes)
  sg2m <- (unique(sort(sg2m)))
  sg2m <- paste("^",sg2m,"$", sep="")
  sg2m <- paste(sg2m,collapse="|")
  Cell.Cycle <- grep(sg2m, x = data$genes, value = TRUE, ignore.case = TRUE)
 # variables
  cellCountLimit = as.numeric(tail(head(data[order(data$meanExp, decreasing = T),],base.mean.rank)[2],1))
  top.rank.line = as.numeric(log2(cellCountLimit))
  SDlimit = dispersion.limit
# add data colors
  data <- data %>%
    mutate(color = ifelse(data$SDs > SDlimit & data$numberOfCells < cellCountLimit,
                          yes = "righLimit",
                          no = ifelse(data$numberOfCells > cellCountLimit,
                                      yes = "leftLimit",
                                      no = "none")))
  #  plot
  myPlot <- ggplot(data,aes(x=log2(numberOfCells),y=SDs, text = genes, label = genes)) +
    geom_point(aes(color = factor(color)),size = cell.size, alpha = cell.transparency) +
    theme_bw(base_size = 16) +
    theme(legend.position = "none") +
    ggtitle(label = "Dispersion Plot", subtitle = "modeled by basemean, number of cells and dispersion") +  # add title
    xlab("log2(number of cells per gene)") +
    ylab("dispersion") +
    geom_vline(xintercept = top.rank.line, colour = rank.line.col) +
    geom_hline(yintercept = SDlimit, colour = disp.line.col) +
    scale_color_manual(values = c("righLimit" = right.sig.col,
                                  "leftLimit" = left.sig.col,
                                  "none" = non.sig.col)) +
    scale_y_continuous(trans = "log1p")
  # geom_text(aes(label=ifelse(genes  %in% mito.genes ,as.character(mito.genes),'')))
  if (mark.mito == TRUE) {
    myPlotGG <- myPlot + geom_text_repel(data = top_labelled,
                                       mapping = aes(label = genes),
                                       size = 3,
                                       fontface = 'bold',
                                       color = 'black',
                                       box.padding = unit(0.5, "lines"),
                                       point.padding = unit(0.5, "lines"))

  }
# get model genes
  my.clust.genes = subset(data, color != "none")[1]
# exclude mito genes from model genes
  if (no.mito.model == TRUE) {
    my.clust.genes <- subset(my.clust.genes, !(genes %in% mito.genes))
  }
  # exclude cell cycle genes from model genes
  if (no.cell.cycle == TRUE) {
    my.clust.genes <- subset(my.clust.genes, !(genes %in% Cell.Cycle))
  }
  ######## if conditions
  #  top_labelled <- subset(data, data$gene %in% Cell.Cycle)
  if (length(Myconds) > 1) {
    for (i in Myconds) {
      data <- subset(DATA, DATA$condition == i)
      data1 <- subset(data, data$SDs > dispersion.limit)
      data1 <- as.character(data1$genes)
      data2 <- head(data[order(data$meanExp, decreasing = TRUE),],base.mean.rank)
      data2 <- as.character(data2$genes)
      bestGenes <- unique(sort(c(data2,data1)))
      if (no.mito.model == TRUE) {
        bestGenes <- subset(bestGenes, !(bestGenes %in% mito.genes))
      }
      # exclude cell cycle genes from model genes
      if (no.cell.cycle == TRUE) {
        bestGenes <- subset(bestGenes, !(bestGenes %in% Cell.Cycle))
      }
      # name
      bestGenes <- as.data.frame(bestGenes)
      NameCol=paste("MyGeneStat",i,sep="_")
      eval(call("<-", as.name(NameCol), bestGenes))
    }
    ## merge
    filenames <- ls(pattern="MyGeneStat_")
    datalist <- mget(filenames)
    Table <- do.call(rbind.data.frame, datalist)
    Table <- as.data.frame(table(as.character(Table$bestGenes)))
    Table <- subset(Table, Table$Freq == length(Myconds))
    my.clust.genes.conds <- as.character(Table$Var1)
  }
  ############# end of loop
  if (length(Myconds) > 1) {
    my.clust.genes <- my.clust.genes.conds
  }
  # retrun
  # write out gene model
#  write.table((my.clust.genes),file="my_model_genes.txt", row.names =FALSE, quote = FALSE, col.names = FALSE)
  gene.counts <- dim(my.clust.genes)[1]
#  message("my_model_genes.txt file is generated, which can be used for clustering.")
if (my.out.put == "plot") {
if (interactive == TRUE) {
  OUT.PUT <- paste(out.name, ".html", sep="")
  htmlwidgets::saveWidget(ggplotly(myPlot),OUT.PUT)
}
if (interactive == FALSE) {
  return(myPlotGG)
 }
}
  if (my.out.put == "data") {
    data <- as.character(as.matrix(my.clust.genes))
    attributes(x)$gene.model <- head(data, gene.num.max)
    return(x)
  }
}
