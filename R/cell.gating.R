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
#' @import shiny
#' @export
cell.gating <- function(x = NULL,
                         my.plot = NULL,
                         plot.type = "tsne") {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  if ("gg" != class(my.plot)[1]) {
    stop("x should be a scatter plot object")
  }
  # geth the genes and scale them based on model
  if (plot.type == "tsne") {
    MyClusts <- x@tsne.data
  }
  if (plot.type == "umap") {
    MyClusts <- x@umap.data
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
        write.table(gr(),file, row.names =F, quote = F, col.names = F)
      })
  }
  #########################################
#  shinyApp(ui, server)
  return(shinyApp(ui, server))
}
