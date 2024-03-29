% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F0002.R
\name{data.aggregation}
\alias{data.aggregation}
\title{Merge multiple data frames and add the condition names to their cell ids}
\usage{
data.aggregation(samples = NULL, condition.names = NULL)
}
\arguments{
\item{samples}{A character vector of data.frame object names.}

\item{condition.names}{A character vector of data.frame condition names.}
}
\value{
An object of class iCellR
}
\description{
This function takes data frame and merges them while also adding condition names to cell ids..
}
\examples{
demo <- read.table(
        file = system.file('extdata', 'demo_data.txt', package = 'iCellR'),
        as.is = TRUE)

# Lets divide your sample in to 3 samples as if you have 3 samples and want to merge them.
sample1 <- demo[1:30]
sample2 <- demo[31:60]
sample3 <- demo[61:90]

# merge all 3 data and add condition names
demo <- data.aggregation(samples =
        c("sample1","sample2","sample3"),
        condition.names = c("WT","ctrl","KO"))
head(demo)[1:4]

# make iCellR object
myDemo.obj <- make.obj(demo)
}
