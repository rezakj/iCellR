#' Filter cells
#'
#' This function takes an object of class iCellR and filters the raw data based on the number of UMIs, genes per cell, percentage of mitochondrial genes per cell, genes, gene expression and cell ids.
#' @param x An object of class iCellR.
#' @param min.mito Min rate for mitochondrial gene expression per cell, default = 0.
#' @param max.mito Max rate for mitochondrial gene expression per cell, default = 1.
#' @param min.genes Min number genes per cell, default = 0.
#' @param max.genes Max number genes per cell, default = Inf.
#' @param min.umis Min number UMIs per cell, default = 0.
#' @param max.umis Max number UMIs per cell, default = Inf.
#' @param filter.by.cell.id A character vector of cell ids to be filtered out.
#' @param keep.cell.id A character vector of cell ids to keep.
#' @param filter.by.gene A character vector of gene names to be filtered by thier expression. If more then one gene is defined it would be OR not AND.
#' @param filter.by.gene.exp.min Minimum gene expression to be filtered by the genes set in filter.by.gene, default = 1.
#' @return An object of class iCellR.
#' @examples
#' demo.obj <- cell.filter(demo.obj,
#'                        min.mito = 0,
#'                        max.mito = 0.05 ,
#'                        min.genes = 100,
#'                        max.genes = 2500,
#'                        min.umis = 0,
#'                        max.umis = Inf)
#'
#' message(demo.obj@my.filters)
#'
#' @export
cell.filter <- function (x = NULL,
                         min.mito = 0,
                         max.mito = 1,
                         min.genes = 0,
                         max.genes = Inf,
                         min.umis = 0,
                         max.umis = Inf,
                         filter.by.cell.id = "character",
                         keep.cell.id = "character",
                         filter.by.gene = "character",
                         filter.by.gene.exp.min = 1) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
# take input from
  DATA <- x@raw.data
# filter by gene
  if (filter.by.gene[1] != "character"){
    genExist <- dim(subset(DATA, row.names(DATA) %in% filter.by.gene))[1]
    numberOfgenes <- length(filter.by.gene)
    if (genExist != numberOfgenes) {
      stop("your data lacks the gene/genes in filter.by.gene.")
    }
  gene.filt <- t(subset(DATA, row.names(DATA) %in% filter.by.gene) >= filter.by.gene.exp.min)
#  gene.filted <- as.data.frame(rowSums(gene.filt) > 0L)
  gene.filted <- as.data.frame(apply(gene.filt, 1, any))
  colnames(gene.filted) <- "gene.filter"
  gene.filt.cells <- row.names(subset(gene.filted, gene.filter == F))
  goneCells <- length(gene.filt.cells)
  allCells <- length(colnames(DATA))
  }
  #
  GenesForFilter = paste(filter.by.gene, collapse=",")
  if (filter.by.gene[1] == "character") {
    GenesForFilter = "no cells were filtered based on expression."
  }
# filter by cell id
#  filter.by.cell.id = c("WT_AAACATACAACCAC.1","WT_AAACATTGATCAGC.1")
  ######## filter out ids
  if (filter.by.cell.id[1] != "character"){
    if (keep.cell.id[1] != "character"){
      stop("if filter.by.cell.id is set you can't set keep.cell.id")
    }
    z = filter.by.cell.id %in% colnames(DATA)
    Missing.cell.ids <- filter.by.cell.id[-which(z, arr.ind = T, useNames = T)]
    Missing.cell.ids <- paste(Missing.cell.ids, collapse=",")
    if (length(filter.by.cell.id) != length(z[z==TRUE])) {
      stop(paste("Your dataset lacks the following ids:   ",Missing.cell.ids))
    }
    DATA <- DATA[ , -which(names(DATA) %in% filter.by.cell.id)]
    FiltCellIds <- paste(filter.by.cell.id , collapse=",")
    message(paste("The following cells were filtered out:   ", FiltCellIds))
  }
  if (filter.by.cell.id[1] == "character") {
    FiltCellIds = "no cell ids were filtered."
  }
  ######## keep ids
  if (keep.cell.id[1] != "character"){
    if (filter.by.cell.id[1] != "character"){
      stop("if keep.cell.id is set you can't set filter.by.cell.id")
    }
    z = keep.cell.id %in% colnames(DATA)
    Missing.cell.ids <- filter.by.cell.id[-which(z, arr.ind = T, useNames = T)]
    Missing.cell.ids <- paste(Missing.cell.ids, collapse=",")
    if (length(keep.cell.id) != length(z[z==TRUE])) {
      stop(paste("Your dataset lacks the following ids:   ",Missing.cell.ids))
    }
    DATA <- DATA[ , which(names(DATA) %in% keep.cell.id)]
    FiltCellIds <- paste(keep.cell.id, collapse=",")
    message(paste("The following cells are kept:   ", FiltCellIds))
  }
  if (keep.cell.id[1] == "character") {
    FiltCellIds = "no cell ids to keep."
  }
# filter by mito
  MAXmito <- as.character(subset(x@stats, x@stats$mito.percent > max.mito)$CellIds)
  MINmito <- as.character(subset(x@stats, x@stats$mito.percent < min.mito)$CellIds)
# filter by nGene
  MAXbad <- as.character(subset(x@stats, x@stats$nGenes > max.genes)$CellIds)
  MINbad <- as.character(subset(x@stats, x@stats$nGenes < min.genes)$CellIds)
# filter by UMI
  MAXbadUMI <- as.character(subset(x@stats, x@stats$UMIs > max.umis)$CellIds)
  MINbadUMI <- as.character(subset(x@stats, x@stats$UMIs < min.umis)$CellIds)
# take all bad genes
  BadMito <- c(MAXmito,MINmito)
  BadGenes <- c(MAXbad,MINbad)
  BadUMIs <- c(MAXbadUMI,MINbadUMI)
  BADgenes <- c(BadMito,BadGenes,BadUMIs)
# filter
  if (length(BadMito) == 0) {
    message("No mito filter")
  } else {
#    DATA <- DATA[ , -which(names(DATA) %in% BadMito)]
    message(paste("cells with min mito ratio of",min.mito,"and max mito ratio of",max.mito,"were filtered."))
  }
  if (length(BadGenes) == 0) {
    message("No gene number filter")
  } else {
#    DATA <- DATA[ , -which(names(DATA) %in% BadGenes)]
    message(paste("cells with min genes of",min.genes,"and max genes of",max.genes,"were filtered."))
  }
  if (length(BadUMIs) == 0) {
    message("No UMI number filter")
  } else {
#    DATA <- DATA[ , -which(names(DATA) %in% BadUMIs)]
    message(paste("cells with min UMIs of",min.umis,"and max UMIs of",max.umis,"were filtered."))
  }
####
  if (length(BADgenes) == 0) {
    message("No UMI, mito and gene number filter")
  } else {
    DATA <- DATA[ , -which(names(DATA) %in% BADgenes)]
      }
####
  if (!exists("gene.filt.cells")) {
    message("No cell filter by provided gene/genes")
  } else {
    DATA <- DATA[ , -which(names(DATA) %in% gene.filt.cells)]
    GenesForFilter = paste(filter.by.gene, collapse=",")
    message(paste(goneCells,"cells out of original",allCells,
                " cells were filtered out as their expression was less than",
                filter.by.gene.exp.min,"for",GenesForFilter))
  }
  if (filter.by.cell.id[1] == "character") {
    message("No cell id filter")
  }
# make a filter parameters file
  FilterFile <- paste("Filters are:
  - min mito",min.mito,"max mito",max.mito,"
  - min # of genes",min.genes,"max # of genes",max.genes,"
  - min UMIs",min.umis,"max UMIs",max.umis,"
  - Cells were also filtered by the expression of the following genes (user input):
  ",GenesForFilter,"
  - In addition the following cell ids were set to keep/filter (user input):
  ",FiltCellIds)
# write filter parameters
#  write.table((FilterFile),file="filters_set.txt", row.names =F, quote = F, col.names = F)
#  message("filters_set.txt file has beed generated and includes the filters set for this experiment.")
  # return data
  row.names(DATA) <- gsub("-",".", row.names(DATA))
  attributes(x)$main.data <- DATA
  attributes(x)$my.filters <- FilterFile
  return(x)
}
