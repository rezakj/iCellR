#' Create MA and Volcano plots.
#'
#' This function takes the result of differential expression (DE) analysis and provides MA and volcano plots.
#' @param x A data frame containing differential expression (DE) analysis results.
#' @param sig.value Choose from "pval" or "padj", default = "padj".
#' @param sig.line A number to draw the line for the significant genes based on sig.value type, default = 0.1.
#' @param plot.type Choose from "ma" or "volcano", default = "volcano".
#' @param x.limit A number to set a limit for the x axis.
#' @param y.limit A number to set a limit for the y axis.
#' @param limit.force If set to TRUE the x.limit and y.limit will be forced, default = FALSE.
#' @param scale.ax If set to TRUE the y axis will be scaled to include all the points, default = TRUE.
#' @param dot.size A number for the size of the points in the plot, default = 1.75.
#' @param dot.transparency Color transparency for points in "scatterplot" and "boxplot", default = 0.5.
#' @param dot.col A set of three colors for the points in the volcano plot, default = c("#E64B35","#3182bd","#636363").
#' @param interactive If set to TRUE an interactive HTML file will be created, default = TRUE.
#' @param out.name If "interactive" is set to TRUE, the output name for HTML, default = "plot".
#' @return Plots
#' @examples
#'
#' diff.res <- run.diff.exp(demo.obj, de.by = "clusters", cond.1 = c(1), cond.2 = c(2))
#'
#' volcano.ma.plot(diff.res,
#'               sig.value = "pval",
#'               sig.line = 0.05,
#'               plot.type = "volcano",
#'               interactive = FALSE)
#'
#' volcano.ma.plot(diff.res,
#'              sig.value = "pval",
#'              sig.line = 0.05,
#'              plot.type = "ma",
#'              interactive = FALSE)
#'
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @importFrom methods new
#' @importFrom stats aggregate as.dendrogram cor cor.test dist hclust p.adjust prcomp quantile sd t.test
#' @importFrom utils capture.output packageVersion read.table write.table
#' @importFrom graphics legend par plot
#' @importFrom plotly layout
#' @importFrom ggplot2 ggplot scale_color_discrete scale_colour_gradient scale_fill_gradient2 scale_x_continuous scale_y_continuous scale_y_discrete stat_summary coord_polar element_rect element_text element_blank facet_wrap scale_color_manual geom_hline geom_jitter geom_vline ylab xlab ggtitle theme_bw aes theme geom_bar geom_point geom_boxplot geom_line
#' @export
volcano.ma.plot <- function (x = NULL,
                          sig.value = "padj",
                          sig.line = 0.1,
                          plot.type = "volcano",
                          x.limit = 2,
                          y.limit = 2,
                          limit.force = FALSE,
                          scale.ax = TRUE,
                          dot.size = 1.75,
                          dot.transparency = 0.5,
                          dot.col = c("#E64B35","#3182bd","#636363"),
                          interactive = TRUE,
                          out.name = "plot") {
  # load data
  DATA <- as.data.frame(x)
  # volcano plot
  if (plot.type == "volcano") {
  # sig line and limits
  Pval = -log10(sig.line)
  YaxiS = y.limit
  XaxiS = x.limit
  # make data ready
  if (sig.value == "padj") {
    results <- subset(DATA,select=c(log2FoldChange,padj))
    data <- data.frame(gene = row.names(results),
                       pvalue = -log10(results$padj),
                       lfc = results$log2FoldChange)
    Myylab <- "-log10 (adjusted p-value)"
  }
  if (sig.value == "pval") {
    results <- subset(DATA,select=c(log2FoldChange,pval))
    data <- data.frame(gene = row.names(results),
                       pvalue = -log10(results$pval),
                       lfc = results$log2FoldChange)
    Myylab <- "-log10 (p-value)"
  }
  # Modify dataset to add new coloumn of colors
  data <- data %>%
    mutate(color = ifelse(data$lfc > 0 & data$pvalue > Pval,
                          yes = "Treated",
                          no = ifelse(data$lfc < 0 & data$pvalue > Pval,
                                      yes = "Untreated",
                                      no = "none")))
  # Color corresponds to fold change directionality
  myPLOT <- ggplot(data, aes(x = lfc, y = pvalue, text = gene)) +
    geom_point(aes(color = factor(color)), size = dot.size, alpha = dot.transparency, na.rm = TRUE) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position = "none") + # remove legend
    ggtitle(label = "Volcano Plot", subtitle = "Colored by directionality") +  # add title
    xlab("log2 (Fold Change)") + # x-axis label
    ylab(Myylab) + # y-axis label
    geom_vline(xintercept = 0, colour = "black") + # add line at 0
    geom_hline(yintercept = Pval, colour = "black") + # p(0.05) = 1.3
    scale_color_manual(values = c("Treated" = dot.col[1],
                                  "Untreated" = dot.col[2],
                                  "none" = dot.col[3]))

  # scale
  ww = c(1:100)
  zz = log1p(1:100)
  # finish plot
  if (limit.force == TRUE) {
    myPLOT <- myPLOT + scale_x_continuous(limits = c(-XaxiS, XaxiS))
  }
  if (scale.ax == TRUE){
    myPLOT <- myPLOT + scale_y_continuous(trans = "log1p")
  }
  if (scale.ax == FALSE){
    myPLOT <- myPLOT + scale_y_continuous(limits = c(YaxiS))
  }
}
#### MA plot
  if (plot.type == "ma") {
  # sig line and limits
  Pval = -log10(sig.line)
  YaxiS = y.limit
  XaxiS = x.limit
  # make data ready
  if (sig.value == "padj") {
    results <- subset(DATA,select=c(baseMean,log2FoldChange,padj))
    data <- data.frame(gene = row.names(results),
                       baseMean = results$baseMean,
                       pvalue = -log10(results$padj),
                       lfc = results$log2FoldChange)
  }
  if (sig.value == "pval") {
    results <- subset(DATA,select=c(baseMean,log2FoldChange,pval))
    data <- data.frame(gene = row.names(results),
                       baseMean = results$baseMean,
                       pvalue = -log10(results$pval),
                       lfc = results$log2FoldChange)
  }
###
  data <- data %>%
    mutate(color = ifelse(data$lfc > 0 & data$pvalue > Pval,
                          yes = "Treated",
                          no = ifelse(data$lfc < 0 & data$pvalue > Pval,
                                      yes = "Untreated",
                                      no = "none")))
  # Color corresponds to fold change directionality
  myPLOT <- ggplot(data, aes(x = baseMean, y = lfc, text = gene)) +
    geom_point(aes(color = factor(color)), size = dot.size, alpha = dot.transparency, na.rm = TRUE) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position = "none") + # remove legend
    ggtitle(label = "MA Plot", subtitle = "Colored by directionality") +  # add title
    xlab("Base Mean") + # x-axis label
    ylab("log2 (Fold Change)") + # y-axis label
    geom_vline(xintercept = 0, colour = "black") + # add line at 0
    geom_hline(yintercept = 0, colour = "black") + # p(0.05) = 1.3
    scale_color_manual(values = c("Treated" = dot.col[1],
                                  "Untreated" = dot.col[2],
                                  "none" = dot.col[3]))
  # scale
  ww = c(1:100)
  zz = log1p(1:100)
  # finish plot
  if (limit.force == TRUE) {
    myPLOT <- myPLOT + scale_y_continuous(limits = c(-YaxiS, YaxiS))
  }
  if (scale.ax == TRUE){
    myPLOT <- myPLOT + scale_x_continuous(trans = "log1p")
  }
  if (scale.ax == FALSE){
    myPLOT <- myPLOT + scale_x_continuous(limits = c(XaxiS))
  }
  }
  # return
  if (interactive == TRUE) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}
