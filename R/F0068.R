#' Read 10X image data
#'
#' This function takes 10X image data files and converts them to proper file format for iCellR.
#' @param dir.10x A directory that includes the 10X image files (scalefactors_json.json, tissue_lowres_image.png and tissue_positions_list.csv).
#' @rdname capture.image.10x
#' @return A list object
#' @importFrom png readPNG
#' @importFrom jsonlite fromJSON
#' @importFrom utils read.csv
#' @export
capture.image.10x <- function (dir.10x = NULL) {
  if (!dir.exists(dir.10x)) {
    stop("Directory is not provided. Please provide a standard 10x image directory
         that includes; scalefactors_json.json, tissue_lowres_image.png and tissue_positions_list.csv files.")
  }
###### read data
  My.image <- readPNG(source = file.path(dir.10x,
                                         "tissue_lowres_image.png"))
  MY.scale.factors <- fromJSON(txt = file.path(dir.10x,
                                               "scalefactors_json.json"))
  My.tissue.positions <- read.csv(file = file.path(dir.10x,
                                                "tissue_positions_list.csv"),
                                  header = FALSE, as.is = TRUE, row.names = 1)
  rownames(My.tissue.positions) <- gsub("-",".",rownames(My.tissue.positions))
#######
  data <- list(My.image,My.tissue.positions,MY.scale.factors)
  names(data) <- c("imgae","position","scale")
  return(data)
}

