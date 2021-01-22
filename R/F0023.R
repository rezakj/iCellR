#' Cell gating
#'
#' This function takes an object of class iCellR and a 2D tSNE or UMAP plot and gates around cells to get their ids.
#' @param x An object of class iCellR.
#' @param my.plot The plot to use for gating. Must be a 2D plot.
#' @param plot.type Choose from knetl, umap and tsne, default = NULL.
#' @return An object of class iCellR.
#' @import shiny
#' @importFrom plotly ggplotly layout plot_ly
#' @export
cell.gating <- function(x = NULL,
                         my.plot = NULL,
                         plot.type = "tsne") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  if (is.null(my.plot)) {
    stop("Please provide a plot")
  }
  if ("gg" != class(my.plot)[1]) {
    stop("plot should be a scatter plot object of class gg")
  }
  # geth the genes and scale them based on model
  if (plot.type == "tsne") {
    MyClusts <- x@tsne.data
  }
  if (plot.type == "knetl") {
    MyClusts <- x@knetl.data
  }
  if (plot.type == "umap") {
    MyClusts <- x@umap.data
  }
  if (plot.type != "umap" && plot.type != "tsne" && plot.type != "knetl") {
    stop("plot type should be chosen. Choose from knetl, tsne or umap.")
  }
  data1 <- paste(MyClusts$V1,MyClusts$V2, sep="_")
  MyClusts <- as.data.frame(cbind(data1,row.names(MyClusts)))
  MyPlot <- ggplotly(my.plot) %>% layout(dragmode = "lasso")
  ##### Make fluid page
  ui <- fluidPage(
    plotlyOutput("plot"),
    downloadButton("downloadData"),
    verbatimTextOutput("click"),
    verbatimTextOutput("brush"))
  #### make server
  server <- function(input, output, session) {
    output$plot <- renderPlotly({
      MyPlot
    })
    #### selcted output
    output$brush <- renderPrint({
      d <- event_data("plotly_selected")
      frame2 <- d
      data1 <- paste(frame2$x,frame2$y, sep="_")
      frame2 <- as.data.frame(data1)
      frame2 <- merge(frame2, MyClusts, by="data1")[2]
      frame2 <- as.character(as.matrix(frame2))
      frame2
    })
    ####### downlaod
    gr <- reactive({
      d <- event_data("plotly_selected")
      frame2 <- d
      data1 <- paste(frame2$x,frame2$y, sep="_")
      frame2 <- as.data.frame(data1)
      frame2 <- merge(frame2, MyClusts, by="data1")[2]
      frame2 <- as.character(as.matrix(frame2))
      frame2
    })
    output$downloadData <- downloadHandler(
      filename = "cellGating.txt",
      content = function(file) {
        write.table(gr(),file, row.names =FALSE, quote = FALSE, col.names = FALSE)
      })
  }
  #########################################
#  shinyApp(ui, server)
  return(shinyApp(ui, server))
}
