% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F0001.R
\name{load10x}
\alias{load10x}
\title{Load 10X data as data.frame}
\usage{
load10x(dir.10x = NULL, gene.name = 2)
}
\arguments{
\item{dir.10x}{A directory that includes the 10X barcodes.tsv, genes.tsv and matrix.mtx files.}

\item{gene.name}{Gene names or ids column number, default = 2.}
}
\value{
The data frame object
}
\description{
This function takes 10X data files barcodes.tsv, genes.tsv and matrix.mtx and converts them to proper matrix file for iCellR.
}
\examples{
my.data <- load10x(system.file("extdata", "filtered_gene_bc_matrices", package = "iCellR"))

# See first few rows and columns
head(my.data)[1:5]

}
