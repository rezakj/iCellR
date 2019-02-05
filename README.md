# iCellR
iCellR is an interactive R package to works with high-throughput single cell sequencing technologies (i.e scRNA-seq, scVDJ-seq and CITE-seq).

We hope to have an official release with stable functions and complete documentation soon!

For citation please use this link (our manuscript is in preparation): https://github.com/rezakj/iCellR

### Single Cell Sequencing R package (iCellR)

<p align="center">
  <img src="https://github.com/rezakj/iCellR/blob/dev/doc/first.gif" width="400"/>
  <img src="https://github.com/rezakj/iCellR/blob/dev/doc/out2.gif" width="400"/>
  <img src="https://github.com/rezakj/iCellR/blob/dev/doc/out3.gif" width="400"/>
  <img src="https://github.com/rezakj/iCellR/blob/dev/doc/out4.gif" width="400"/> 
<img src="https://github.com/rezakj/iCellR/blob/master/doc/out10.gif" /> 
</p>

***
## How to install iCellR
        
```r
library(devtools)
install_github("rezakj/iCellR")

# or
git clone https://github.com/rezakj/iCellR.git
R
install.packages('iCellR/', repos = NULL, type="source")
```

## Download a sample data

- Download and unzip a publicly available sample [PBMC](https://en.wikipedia.org/wiki/Peripheral_blood_mononuclear_cell) scRNA-Seq data.

```r
# set your working directory 
setwd("/your/download/directory")

# save the URL as an object
sample.file.url = "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"

# download the file
download.file(url = sample.file.url, 
     destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz", 
     method = "auto")  

# unzip the file. 
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")    
```

***
# How to use iCellR for analyzing scRNA-seq data

To run a test sample follow these steps:

- Go to the R environment load the iCellR package and the PBMC sample data that you downloaded.

```r
library("iCellR")
my.data <- load10x("filtered_gene_bc_matrices/hg19/",gene.name = "geneSymbol")
```

To see the help page for each function use question mark as: 

```r
?load10x
```

- Aggregate data
     
Conditions in iCellR are set in the header of the data and are separated by an underscore (_). Let's say you want to merge multiple datasets and run iCellR in aggregate mode. Here’s an example: I divided this sample into three sets and then aggregate them into one matrix. 

```r
dim(my.data)
# [1] 32738  2700

# divide your sample into three samples for this example 
  sample1 <- my.data[1:900]
  sample2 <- my.data[901:1800]
  sample3 <- my.data[1801:2700]
  
# merge all of your samples to make a single aggregated file.    
my.data <- data.aggregation(samples = c("sample1","sample2","sample3"), 
	condition.names = c("WT","KO","Ctrl"))
```

- Check the head of your file.

```r
# here is how the head of the first 2 cells in the aggregated file looks like.	
head(my.data)[1:2]
#         WT_AAACATACAACCAC-1 WT_AAACATTGAGCTAC-1
#A1BG                       0                   0
#A1BG.AS1                   0                   0
#A1CF                       0                   0
#A2M                        0                   0
#A2M.AS1                    0                   0

# as you see the header has the conditions now
```


- Make an object of class iCellR.

```r
my.obj <- make.obj(my.data)
my.obj
#[1] "An object of class iCellR version: 0.99.0"                                     
#[2] "Raw/original data dimentions (rows,columns): 32738,2700"                       
#[3] "Data conditions in raw data: Ctrl,KO,WT (900,900,900)"                         
#[4] "Columns names: WT_AAACATACAACCAC.1,WT_AAACATTGAGCTAC.1,WT_AAACATTGATCAGC.1 ..."
#[5] "Row names: A1BG,A1BG.AS1,A1CF ..."   
```

- Perform some QC 

```r
my.obj <- qc.stats(my.obj)
``` 

- Plot QC

By default all the plotting functions would create interactive html files unless you set this parameter: interactive = FALSE.

```r
# plot UMIs, genes and percent mito all at once and in one plot. 
# you can make them individually as well, see the arguments ?stats.plot.
stats.plot(my.obj,
	plot.type = "all.in.one",
	out.name = "UMI-plot",
	interactive = FALSE,
	cell.color = "slategray3", 
	cell.size = 1, 
	cell.transparency = 0.5,
	box.color = "red",
	box.line.col = "green")
```

<p align="center">
  <img src="https://github.com/rezakj/iCellR/blob/dev/doc/stats.png" />
</p>

```r  
# Scatter plots
stats.plot(my.obj, plot.type = "point.mito.umi", out.name = "mito-umi-plot")
stats.plot(my.obj, plot.type = "point.gene.umi", out.name = "gene-umi-plot")
```
<p align="center">
  <img src="https://github.com/rezakj/iCellR/blob/dev/doc/out5.gif" width="400"/>
  <img src="https://github.com/rezakj/iCellR/blob/dev/doc/out6.gif" width="400"/>
</p>

- Filter cells. 

iCellR allows you to filter based on library sizes (UMIs), number of genes per cell, percent mitochondrial content, one or more genes, and cell ids.

```r
my.obj <- cell.filter(my.obj,
	min.mito = 0,
	max.mito = 0.05,
	min.genes = 200,
	max.genes = 2400,
	min.umis = 0,
	max.umis = Inf)
	
#[1] "cells with min mito ratio of 0 and max mito ratio of 0.05 were filtered."
#[1] "cells with min genes of 200 and max genes of 2400 were filtered."
#[1] "No UMI number filter"
#[1] "No cell filter by provided gene/genes"
#[1] "No cell id filter"
#[1] "filters_set.txt file has beed generated and includes the filters set for this experiment."	

# more examples 
# my.obj <- cell.filter(my.obj, filter.by.gene = c("RPL13","RPL10")) # filter our cell having no counts for these genes
# my.obj <- cell.filter(my.obj, filter.by.cell.id = c("WT_AAACATACAACCAC.1")) # filter our cell cell by their cell ids.

# chack to see how many cells are left.  
dim(my.obj@main.data)
#[1] 32738  2637
```
- Down sampling 

This step is optional and is for having the same number of cells for each condition. 

```r
# optional
# my.obj <- down.sample(my.obj)
#[1] "From"
#[1] "Data conditions: Ctrl,KO,WT (877,877,883)"
#[1] "to"
#[1] "Data conditions: Ctrl,KO,WT (877,877,877)"
```

- Normalize data

