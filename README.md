# iCellR
iCellR is an interactive R package to work with high-throughput single cell sequencing technologies (i.e scRNA-seq, scVDJ-seq and CITE-seq).

Link to Comprehensive R Archive Network [(CRAN)](https://cran.r-project.org/web/packages/iCellR/index.html)

Link to manual: [Manual](https://cran.r-project.org/web/packages/iCellR/iCellR.pdf)

Link to a video tutorial for CITE-Seq and scRNA-Seq analysis: [Video](https://vimeo.com/337822487)

For citation please use this link (our manuscript is in preparation): https://CRAN.R-project.org/package=iCellR 

If you are using FlowJo or SeqGeq, they have made plugins for iCellR and other single cell tools: https://www.flowjo.com/exchange/#/ (list of all plugins) and https://www.flowjo.com/exchange/#/plugin/profile?id=34 (iCellR plugin)

Tutorials: [example 1 code](https://genome.med.nyu.edu/results/external/iCellR/example1/code.txt) and [results](https://genome.med.nyu.edu/results/external/iCellR/example1/)

### Single (i) Cell R package (iCellR)

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/first.gif" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/out2.gif" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/out3.gif" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/out4.gif" width="400"/> 
	  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/gating2.gif"/>
</p>

***
## How to install iCellR
        
```r
# Install from CRAN 
install.packages("iCellR")

# Install from github
#library(devtools)
#install_github("rezakj/iCellR")

# or
#git clone https://github.com/rezakj/iCellR.git
#R
#install.packages('iCellR/', repos = NULL, type="source")
```

## Download a sample data

- Download and unzip a publicly available sample [PBMC](https://en.wikipedia.org/wiki/Peripheral_blood_mononuclear_cell) scRNA-Seq data.

```r
# set your working directory 
setwd("/your/download/directory")

# save the URL as an object
sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/pbmc3k_filtered_gene_bc_matrices.tar.gz"

# download the file
download.file(url = sample.file.url, 
     destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz", 
     method = "auto")  

# unzip the file. 
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")
```
more data available here: 
https://genome.med.nyu.edu/results/external/iCellR/

***
# How to use iCellR for analyzing scRNA-seq data

To run a test sample follow these steps:

- Go to the R environment load the iCellR package and the PBMC sample data that you downloaded.

```r
library("iCellR")
my.data <- load10x("filtered_gene_bc_matrices/hg19/")

# This directory includes; barcodes.tsv, genes.tsv/features.tsv and matrix.mtx files (data could be zipped or unzipped).

# if your data is in a csv or tsv format read it like this example
# my.data <- read.delim("CITE-Seq_sample_RNA.tsv.gz",header=TRUE)
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
###################################
,--. ,-----.       ,--.,--.,------.
`--''  .--./ ,---. |  ||  ||  .--. '
,--.|  |    | .-. :|  ||  ||  '--'.'
|  |'  '--'\   --. |  ||  ||  |
`--' `-----' `----'`--'`--'`--' '--'
###################################
An object of class iCellR version: 0.99.0 
    Raw/original data dimentions (rows,columns): 32738,2700 
    Data conditions in raw data: Ctrl,KO,WT (900,900,900) 
    Row names: A1BG,A1BG.AS1,A1CF ... 
    Columns names: WT_AAACATACAACCAC.1,WT_AAACATTGAGCTAC.1,WT_AAACATTGATCAGC.1 ... 
###################################
   QC stats performed: FALSE , PCA performed: FALSE , CCA performed: FALSE
   Clustering performed: FALSE , Number of clusters: 0
   tSNE performed: FALSE , UMAP performed: FALSE , DiffMap performed: FALSE
   Main data dimentions (rows,columns): 0 0
   Normalization factors:  ...
   Imputed data dimentions (rows,columns): 0 0
############## scVDJ-Seq ###########
   VDJ data dimentions (rows,columns): 0 0
############## CITE-Seq ############
   ADT raw data dimentions (rows,columns): 0 0
   ADT main data dimentions (rows,columns): 0 0
   ADT columns names:  ...
   ADT row names:  ...
######## iCellR object made ########
```

- Perform some QC 

```r
my.obj <- qc.stats(my.obj)
``` 


- Cell cycle prediction 

```r
my.obj <- cc(my.obj, s.genes = s.phase, g2m.genes = g2m.phase)
head(my.obj@stats)

#                                CellIds nGenes UMIs mito.percent
#WT_AAACATACAACCAC.1 WT_AAACATACAACCAC.1    781 2421  0.030152829
#WT_AAACATTGAGCTAC.1 WT_AAACATTGAGCTAC.1   1352 4903  0.037935958
#WT_AAACATTGATCAGC.1 WT_AAACATTGATCAGC.1   1131 3149  0.008891712
#WT_AAACCGTGCTTCCG.1 WT_AAACCGTGCTTCCG.1    960 2639  0.017430845
#WT_AAACCGTGTATGCG.1 WT_AAACCGTGTATGCG.1    522  981  0.012232416
#WT_AAACGCACTGGTAC.1 WT_AAACGCACTGGTAC.1    782 2164  0.016635860
#                    S.phase.probability g2m.phase.probability      S.Score
#WT_AAACATACAACCAC.1        0.0012391574          0.0004130525  0.030569081
#WT_AAACATTGAGCTAC.1        0.0002039568          0.0004079135 -0.077860621
#WT_AAACATTGATCAGC.1        0.0003175611          0.0019053668 -0.028560560
#WT_AAACCGTGCTTCCG.1        0.0007578628          0.0011367942  0.001917225
#WT_AAACCGTGTATGCG.1        0.0000000000          0.0020387360 -0.020085210
#WT_AAACGCACTGGTAC.1        0.0000000000          0.0000000000 -0.038953135
#                        G2M.Score Phase
#WT_AAACATACAACCAC.1 -0.0652390011     S
#WT_AAACATTGAGCTAC.1 -0.1277015099    G1
#WT_AAACATTGATCAGC.1 -0.0036505733    G1
#WT_AAACCGTGCTTCCG.1 -0.0499511543     S
#WT_AAACCGTGTATGCG.1  0.0009426363   G2M
#WT_AAACGCACTGGTAC.1 -0.0680240629    G1


# plot cell cycle rate
pie(table(my.obj@stats$Phase))
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/iCellR_1.png" width="400"/>
</p>

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
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/stats.png" />
</p>

```r  
# Scatter plots
stats.plot(my.obj, plot.type = "point.mito.umi", out.name = "mito-umi-plot")
stats.plot(my.obj, plot.type = "point.gene.umi", out.name = "gene-umi-plot")
```
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/out5.gif" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/out6.gif" width="400"/>
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

You have a few options to normalize your data based on your study. You can also normalize your data using tools other than iCellR and import your data to iCellR. We recommend "ranked.glsf" normalization for most single cell studies. This normalization is great for fixing matrixes with lots of zeros and because it's geometric it is great for fixing for batch effects, as long as all the data is aggregated into one file (to aggregate your data see "aggregating data" section above). 

```r
my.obj <- norm.data(my.obj, 
     norm.method = "ranked.glsf",
     top.rank = 500) # best for scRNA-Seq

# more examples
#my.obj <- norm.data(my.obj, norm.method = "ranked.deseq", top.rank = 500)
#my.obj <- norm.data(my.obj, norm.method = "deseq") # best for bulk RNA-Seq 
#my.obj <- norm.data(my.obj, norm.method = "global.glsf") # best for bulk RNA-Seq 
#my.obj <- norm.data(my.obj, norm.method = "rpm", rpm.factor = 100000) # best for bulk RNA-Seq
#my.obj <- norm.data(my.obj, norm.method = "spike.in", spike.in.factors = NULL)
#my.obj <- norm.data(my.obj, norm.method = "no.norm") # if the data is already normalized
```

- Perform second QC 

```r
my.obj <- qc.stats(my.obj,which.data = "main.data")

stats.plot(my.obj,
	plot.type = "all.in.one",
	out.name = "UMI-plot",
	interactive = F,
	cell.color = "slategray3", 
	cell.size = 1, 
	cell.transparency = 0.5,
	box.color = "red",
	box.line.col = "green",
	back.col = "white")
``` 

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/stats2.png" />
</p>

- Scale data (optional)

iCellR dose not need this step as it scales the data when they need to be scaled on the fly; like for plotting or PCA. 
It is important to use the untansformed data for differential expression analysis to calculate the accurate fold changes.
If you run this function the scaled data will be saved in different slot for you to download for plotting but will not be use by iCellR.

```r
# my.obj <- data.scale(my.obj)
```

- Gene stats

```r
my.obj <- gene.stats(my.obj, which.data = "main.data")

head(my.obj@gene.data[order(my.obj@gene.data$numberOfCells, decreasing = T),])
#       genes numberOfCells totalNumberOfCells percentOfCells  meanExp
#30303 TMSB4X          2637               2637      100.00000 38.55948
#3633     B2M          2636               2637       99.96208 45.07327
#14403 MALAT1          2636               2637       99.96208 70.95452
#27191 RPL13A          2635               2637       99.92416 32.29009
#27185  RPL10          2632               2637       99.81039 35.43002
#27190  RPL13          2630               2637       99.73455 32.32106
#               SDs condition
#30303 7.545968e-15       all
#3633  2.893940e+01       all
#14403 7.996407e+01       all
#27191 2.783799e+01       all
#27185 2.599067e+01       all
#27190 2.661361e+01       all
```

- Make a gene model for clustering

It's best to always to avoid global clustering and use a set of model genes. In bulk RNA-seq data it is very common to cluster the samples based on top 500 genes ranked by base mean, this is to reduce the noise. In scRNA-seq data, it's great to do so as well. This coupled with our ranked.glsf normalization is good for matrices with a lot of zeros. You can also use your set of genes as a model rather than making one. 

```r
# See model plot 
make.gene.model(my.obj, my.out.put = "plot",
	dispersion.limit = 1.5, 
	base.mean.rank = 500, 
	no.mito.model = T, 
	mark.mito = T, 
	interactive = F,
	no.cell.cycle = T,
	out.name = "gene.model")
	
# Write the gene model data into the object

my.obj <- make.gene.model(my.obj, my.out.put = "data",
	dispersion.limit = 1.5, 
	base.mean.rank = 500, 
	no.mito.model = T, 
	mark.mito = T, 
	interactive = F,
	no.cell.cycle = T,
	out.name = "gene.model")

head(my.obj@gene.model)
# "ACTB"  "ACTG1" "ACTR3" "AES"   "AIF1"  "ALDOA"

# get html plot (optional)
#make.gene.model(my.obj, my.out.put = "plot",
#	dispersion.limit = 1.5, 
#	base.mean.rank = 500, 
#	no.mito.model = T, 
#	mark.mito = T, 
#	interactive = T,
#	out.name = "plot4_gene.model")
```
To view an the html intractive plot click on this links: [Dispersion plot](https://rawgit.com/rezakj/scSeqR/dev/doc/gene.model.html)


<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/gene.model.png" width="800" height="800" />
</p>


- Perform Principal component analysis (PCA)

Note: skip this step if you did batch correction. For batch correction (sample alignment/harmonization) see the sections; CCA, MNN or anchor (integration) alignment. 

```r
my.obj <- run.pca(my.obj, method = "gene.model", gene.list = my.obj@gene.model,data.type = "main")

opt.pcs.plot(my.obj)

# 2 round PCA (to find top genes in the first 10 PCs and re-run PCA for better clustering
## This is optional and might not be good in some cases
length(my.obj@gene.model)
# 681
my.obj <- find.dim.genes(my.obj, dims = 1:10,top.pos = 20, top.neg = 20) # (optional)

length(my.obj@gene.model)
# 214

# second round PC
my.obj <- run.pca(my.obj, method = "gene.model", gene.list = my.obj@gene.model,data.type = "main")
```        

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/Opt_Number_Of_PCs.png" />
</p>


- Cluster the data

Here we cluster the first 10 dimensions of the data which is converted to principal components. You have the option of clustering your data based on the following methods: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans"

 For the distance calculation used for clustering, you have the following options: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski" or "NULL"

 With the following indexing methods: "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw"

We recomand to use the defult options as below:

```r
my.obj <- run.clustering(my.obj, 
	clust.method = "kmeans", 
	dist.method = "euclidean",
	index.method = "silhouette",
	max.clust = 25,
	min.clust = 2,
	dims = 1:10)

# or 

library(Rphenograph)
my.obj <- run.phenograph(my.obj,k = 100,dims = 1:10)


# if Rphenograph not installed
#devtools::install_github("JinmiaoChenLab/Rphenograph")

# If you want to manually set the number of clusters, and not used the predicted optimal number, set the minimum and maximum to the number you want:
#my.obj <- run.clustering(my.obj, 
#	clust.method = "ward.D",
#	dist.method = "euclidean",
#	index.method = "ccc",
#	max.clust = 8,
#	min.clust = 8,
#	dims = 1:10)

# more examples 

#my.obj <- run.clustering(my.obj, 
#	clust.method = "ward.D", 
#	dist.method = "euclidean",
#	index.method = "kl",
#	max.clust = 25,
#	min.clust = 2,
#	dims = 1:10)
```

- Perform Dimensionality reduction

```r
# tSNE
my.obj <- run.pc.tsne(my.obj, dims = 1:10)

# UMAP
my.obj <- run.umap(my.obj, dims = 1:10, method = "naive")
# or 
# my.obj <- run.umap(my.obj, dims = 1:10, method = "umap-learn") 

# diffusion map
# this requires python packge phate 
# pip install --user phate
# Install phateR version 2.9
# wget https://cran.r-project.org/src/contrib/Archive/phateR/phateR_0.2.9.tar.gz
# install.packages('phateR/', repos = NULL, type="source")

# optional 
# library(phateR)
# my.obj <- run.diffusion.map(my.obj, dims = 1:10, method = "phate")
```

- Visualize data

```r
# clusters
A= cluster.plot(my.obj,plot.type = "pca",interactive = F)
B= cluster.plot(my.obj,plot.type = "umap",interactive = F)
C= cluster.plot(my.obj,plot.type = "tsne",interactive = F) 
D= cluster.plot(my.obj,plot.type = "diffusion",interactive = F)

library(gridExtra)
grid.arrange(A,B,C,D)

# conditions 
A= cluster.plot(my.obj,plot.type = "pca",col.by = "conditions",interactive = F)
B= cluster.plot(my.obj,plot.type = "umap",col.by = "conditions",interactive = F)
C= cluster.plot(my.obj,plot.type = "tsne",col.by = "conditions",interactive = F)
D= cluster.plot(my.obj,plot.type = "diffusion",col.by = "conditions",interactive = F)

library(gridExtra)
grid.arrange(A,B,C,D)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/1_AllClusts.png"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/2_AllConds.png"/>      
</p>

- 3D plots, density plots and intractive plots 

```r
# 2D
cluster.plot(my.obj,
	cell.size = 1,
	plot.type = "tsne",
	cell.color = "black",
	back.col = "white",
	col.by = "clusters",
	cell.transparency = 0.5,
	clust.dim = 2,
	interactive = F)
	
# intractive 2D
cluster.plot(my.obj,
	plot.type = "tsne",
	col.by = "clusters",
	clust.dim = 2,
	interactive = T,
	out.name = "tSNE_2D_clusters")

# intractive 3D
cluster.plot(my.obj,
	plot.type = "tsne",
	col.by = "clusters",
	clust.dim = 3,
	interactive = T,
	out.name = "tSNE_3D_clusters")

# Density plot for clusters 
cluster.plot(my.obj,
	plot.type = "pca",
	col.by = "clusters",
	interactive = F,
	density=T)

# Density plot for conditions 
cluster.plot(my.obj,
	plot.type = "pca",
	col.by = "conditions",
	interactive = F,
	density=T)
```
## To see the above made interactive plots click on these links: [2Dplot](https://rawgit.com/rezakj/scSeqR/dev/doc/tSNE_2D_clusters.html) and [3Dplot](https://rawgit.com/rezakj/scSeqR/dev/doc/tSNE_3D_clusters.html)
        
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_2D_clusters.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_3D.png" width="400"/> 
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/density_conditions.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/density_clusters.png" width="400"/> 	
</p>

- More plots

```r
# plot 
cluster.plot(my.obj,
	cell.size = 1,
	plot.type = "diffusion",
	cell.color = "black",
	back.col = "white",
	col.by = "clusters",
	cell.transparency = 0.5,
	clust.dim = 2,
	interactive = F)
	
cluster.plot(my.obj,
	cell.size = 1,
	plot.type = "diffusion",
	cell.color = "black",
	back.col = "white",
	col.by = "clusters",
	cell.transparency = 0.5,
	clust.dim = 3,
	interactive = F)	
	
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/Diffusion.png" width="400"/>
	<img src="https://github.com/rezakj/scSeqR/blob/dev/doc/diffiusion3D.gif" width="400"/>
</p>

- Plotting clusters and conditions at the same time

```r
cluster.plot(my.obj,
	cell.size = 1,
	plot.type = "umap",
	cell.color = "black",
	back.col = "white",
	col.by = "clusters",
	cell.transparency = 0.5,
	clust.dim = 2,
	cond.shape = T,
	interactive = T,
	out.name = "2d_UMAP_clusters_conds")
```


<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/Conds_clusts.gif"/>
</p>


- Normalized cell frequencies in clusters and conditions

```r
# If normalize.ncell = TRUE it would down sample the conditions randomly so all the conditions have equal number of cells, if FALSE it would output the raw cell counts.

# bar plot
clust.cond.info(my.obj, plot.type = "bar", normalize.ncell = FALSE, my.out.put = "plot")
# Pie chart 
clust.cond.info(my.obj, plot.type = "pie", normalize.ncell = FALSE, ,my.out.put = "plot")

# data 
my.obj <- clust.cond.info(my.obj, plot.type = "bar", normalize.ncell = F)
#head(my.obj@my.freq)
#  conditions clusters Freq
#1       ctrl        1  199
#2         KO        1  170
#3         WT        1  182
#4       ctrl        2  106
#5         KO        2  116
#6         WT        2  113
```
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/3_clust_cond_freq_info_bar.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/4_clust_cond_freq_info_pie.png" width="400"/> 	
</p>

- Avrage expression per cluster

```r
# remember to run this command again if you recluster the data (change the number of clusters). 
# This is because this data slot with avrage expressions per cluster is used for finding markers and downstream analysis. 

my.obj <- clust.avg.exp(my.obj)

head(my.obj@clust.avg)
#      gene   cluster_1   cluster_2  cluster_3 cluster_4   cluster_5   cluster_6
#1     A1BG 0.072214723 0.092648973 0.08258609         0 0.027183115 0.072291636
#2 A1BG.AS1 0.014380756 0.003280237 0.01817982         0 0.000000000 0.011545546
#3     A1CF 0.000000000 0.000000000 0.00000000         0 0.000000000 0.000000000
#4      A2M 0.000000000 0.000000000 0.00000000         0 0.007004131 0.004672857
#5  A2M.AS1 0.003520828 0.039985296 0.00876364         0 0.056596203 0.018445562
#6    A2ML1 0.000000000 0.000000000 0.00000000         0 0.000000000 0.000000000
#   cluster_7  cluster_8   cluster_9
#1 0.09058946 0.04466827 0.027927923
#2 0.00000000 0.01534541 0.005930566
#3 0.00000000 0.00000000 0.000000000
#4 0.00000000 0.00000000 0.003411938
#5 0.00000000 0.00000000 0.000000000
#6 0.00000000 0.00000000 0.000000000
```

- Save your object

```r
save(my.obj, file = "my.obj.Robj")
```        

- Find marker genes

```r
marker.genes <- findMarkers(my.obj,
	fold.change = 2,
	padjval = 0.1)

dim(marker.genes)
# [1] 1070   17

head(marker.genes)
#                row   baseMean    baseSD AvExpInCluster AvExpInOtherClusters
#LRRN3         LRRN3 0.01428477 0.1282046     0.05537243          0.003437002
#LINC00176 LINC00176 0.06757573 0.2949763     0.21404151          0.028906516
#FHIT           FHIT 0.10195359 0.3885343     0.31404936          0.045957058
#TSHZ2         TSHZ2 0.04831334 0.2628778     0.14300998          0.023311970
#CCR7           CCR7 0.28132627 0.6847417     0.81386444          0.140728033
#SCGB3A1     SCGB3A1 0.06319598 0.3554273     0.18130557          0.032013232
#          foldChange log2FoldChange         pval         padj clusters
#LRRN3      16.110677       4.009945 1.707232e-06 2.847662e-03        1
#LINC00176   7.404611       2.888424 4.189197e-16 7.117446e-13        1
#FHIT        6.833539       2.772633 1.576339e-19 2.681353e-16        1
#TSHZ2       6.134616       2.616973 8.613622e-10 1.455702e-06        1
#CCR7        5.783243       2.531879 1.994533e-42 3.400679e-39        1
#SCGB3A1     5.663457       2.501683 2.578484e-07 4.313805e-04        1
#               gene  cluster_1   cluster_2   cluster_3 cluster_4   cluster_5
#LRRN3         LRRN3 0.05537243 0.004102916 0.002190847         0 0.010902326
#LINC00176 LINC00176 0.21404151 0.016772401 0.005203161         0 0.009293024
#FHIT           FHIT 0.31404936 0.008713243 0.022934924         0 0.035701186
#TSHZ2         TSHZ2 0.14300998 0.008996236 0.009444180         0 0.000000000
#CCR7           CCR7 0.81386444 0.075719109 0.034017494         0 0.021492756
#SCGB3A1     SCGB3A1 0.18130557 0.039644151 0.001183264         0 0.000000000
#            cluster_6  cluster_7   cluster_8   cluster_9
#LRRN3     0.002087831 0.00000000 0.000000000 0.012113258
#LINC00176 0.086762509 0.01198777 0.003501552 0.003560614
#FHIT      0.104189143 0.04144293 0.041064681 0.007218861
#TSHZ2     0.065509372 0.01690584 0.002352707 0.015350123
#CCR7      0.272580821 0.06523324 0.257130255 0.031304151
#SCGB3A1   0.078878071 0.01198777 0.000000000 0.043410608

# baseMean: average expression in all the cells
# baseSD: Standard Deviation
# AvExpInCluster: average expression in cluster number (see clusters)
# AvExpInOtherClusters: average expression in all the other clusters
# foldChange: AvExpInCluster/AvExpInOtherClusters
# log2FoldChange: log2(AvExpInCluster/AvExpInOtherClusters)
# pval: P value 
# padj: Adjusted P value 
# clusters: marker for cluster number
# gene: marker gene for the cluster
# the rest are the average expression for each cluster
```

- Plot genes

```r
# tSNE 2D
A <- gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
# PCA 2D	
B <- gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "pca")
	
# Box Plot
C <- gene.plot(my.obj, gene = "MS4A1", 
	box.to.test = 0, 
	box.pval = "sig.signs",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
D <- gene.plot(my.obj, gene = "MS4A1", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
	
library(gridExtra)
grid.arrange(A,B,C,D)	
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/7_genePlots.png"/>
</p>


- Multiple plots

```r
genelist = c("PPBP","LYZ","MS4A1","GNLY","LTB","NKG7","IFITM2","CD14","S100A9")
###
for(i in genelist){
	MyPlot <- gene.plot(my.obj, gene = i, 
		interactive = F,
		plot.data.type = "umap",
		cell.transparency = 1)
	NameCol=paste("PL",i,sep="_")
	eval(call("<-", as.name(NameCol), MyPlot))
}
###
library(cowplot)
filenames <- ls(pattern="PL_")
plot_grid(plotlist=mget(filenames[1:9]))
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/8_genePlots.png" />
</p>

- Heatmap

```r
# find top genes
MyGenes <- top.markers(marker.genes, topde = 10, min.base.mean = 0.2,filt.ambig = F)
MyGenes <- unique(MyGenes)
# plot
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = T, out.name = "plot", cluster.by = "clusters")
# or 
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters")
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/9_heatmap_gg.png" />
</p>

 - Run data imputation

```r
my.obj <- run.impute(my.obj, dims = 1:10, cell.ratio = 2, data.type = "pca")

# more examples
# my.obj <- run.impute(my.obj, cell.ratio = 2, data.type = "tsne")
# my.obj <- run.impute(my.obj, cell.ratio = 2, data.type = "umap")

# save after imputation 
save(my.obj, file = "my.obj.Robj")

# some more plots from another analysis 
A=heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", cell.sort = TRUE)
B=heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", data.type = "imputed", cell.sort = TRUE)
C=heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "conditions", cell.sort = TRUE)
D=heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "none", data.type = "imputed", cell.sort = TRUE)

# If cluster.by = "none", your heamap would have like a Pseudotime effect.
# It calculates the distance between the cells based on the genes in the heatmap. 

library(gridExtra)
grid.arrange(A,B,C,D)

# main data 
gene.plot(my.obj, gene = "MS4A1", 
    plot.type = "scatterplot",
    interactive = F,
    data.type = "main")

# imputed data 
gene.plot(my.obj, gene = "MS4A1", 
    plot.type = "scatterplot",
    interactive = F,
    data.type = "imputed")		
```

<p align="center">
	<img src="https://github.com/rezakj/scSeqR/blob/master/doc/heatmaps.png" />
	<img src="https://github.com/rezakj/scSeqR/blob/dev/doc/imputed_dotPlot.png" />
	<img src="https://github.com/rezakj/scSeqR/blob/dev/doc/imputed_BoxPlot.png" />
</p>

 - Plotting conditions and clusters for genes
 
 ```r
 A <- gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	cell.transparency = 1,
	scaleValue = TRUE,
	min.scale = 0,
	max.scale = 2.5,
	back.col = "white",
	cond.shape = TRUE)
B <- gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	cell.transparency = 1,
	scaleValue = TRUE,
	min.scale = 0,
	max.scale = 2.5,
	back.col = "white",
	cond.shape = TRUE,
	conds.to.plot = c("KO","WT"))

C <- gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "boxplot",
	interactive = F,
	back.col = "white",
	cond.shape = TRUE,
	conds.to.plot = c("KO"))

D <- gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "barplot",
	interactive = F,
	cell.transparency = 1,
	back.col = "white",
	cond.shape = TRUE,
	conds.to.plot = c("KO","WT"))

library(gridExtra)
grid.arrange(A,B,C,D)
 ```
 
 <p align="center">
	<img src="https://github.com/rezakj/scSeqR/blob/master/doc/Conds.png" />
</p>
 
 
 - Cell type prediction using ImmGen
 
Note that ImmGen is mouse genome data and the sample data here is human. For 157 ULI-RNA-Seq samples use this meta data: [metadata](https://github.com/rezakj/scSeqR/blob/dev/doc/uli_RNA_metadat.txt). 

```r
Cluster = 8
MyGenes <- top.markers(marker.genes, topde = 40, min.base.mean = 0.2, cluster = Cluster)
# plot 
cell.type.pred(immgen.data = "rna", gene = MyGenes, plot.type = "point.plot")
# and
cell.type.pred(immgen.data = "uli.rna", gene = MyGenes, plot.type = "point.plot", top.cell.types = 50)
# or 
cell.type.pred(immgen.data = "rna", gene = MyGenes, plot.type = "heatmap")
# and
cell.type.pred(immgen.data = "uli.rna", gene = MyGenes, plot.type = "heatmap")

# And finally check the genes in the cells and find the common ones to predict
# heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters")

# As you can see cluster 8 is most likely to be B-cells.  

# for tissue type prediction use this:
#cell.type.pred(immgen.data = "mca", gene = MyGenes, plot.type = "point.plot")
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/ImmGen_pointPlot_RNA_Cluster_7.png" />
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/ImmGen_pointPlot_ULI-RNA_Cluster_7.png" /> 
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/ImmGen_heatmap_RNA_Cluster_7.png" /> 
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/ImmGen_heatmap_ULI-RNA_Cluster_7.png" />	
</p>


 - Pathway analysis
 
```r
# Pathway  
# pathways.kegg(my.obj, clust.num = 7) 
# this function is being improved and soon will be available
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/7_cluster_KEGGpathways.png" />    
</p>

- QC on clusters 

```r
clust.stats.plot(my.obj, plot.type = "box.mito", interactive = F)
clust.stats.plot(my.obj, plot.type = "box.gene", interactive = F)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/6_cluster_mito_ratio.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/5_cluster_gene_cov.png" width="400"/>      
</p>


- Differential Expression Analysis 

The differential expression (DE) analysis function in iCellR allows the users to choose from any combinations of clusters and conditions. For example, a user with two samples (say WT and KO) has four different possible ways of comparisons:

a-Comparing a cluster/clusters with different cluster/clusters (e.g. cluster 1 and 2 vs. 4)

b-Comparing a cluster/clusters with different cluster/clusters only in one/more condition/conditions (e.g. cluster 1 vs cluster 2 but only the WT sample)

c-Comparing a condtion/condtions with different condtion/condtions (e.g. WT vs KO)

d-Comparing a condtion/condtions with different condtion/condtions only in one/more cluster/clusters (e.g. cluster 1 WT vs cluster 1 KO)

```r
diff.res <- run.diff.exp(my.obj, de.by = "clusters", cond.1 = c(1,4), cond.2 = c(2))
diff.res1 <- as.data.frame(diff.res)
diff.res1 <- subset(diff.res1, padj < 0.05)
head(diff.res1)
#             baseMean        1_4           2 foldChange log2FoldChange         pval
#AAK1       0.19554589 0.26338228 0.041792762 0.15867719      -2.655833 8.497012e-33
#ABHD14A    0.09645732 0.12708519 0.027038379 0.21275791      -2.232715 1.151865e-11
#ABHD14B    0.19132829 0.23177944 0.099644572 0.42991118      -1.217889 3.163623e-09
#ABLIM1     0.06901900 0.08749258 0.027148089 0.31029018      -1.688310 1.076382e-06
#AC013264.2 0.07383608 0.10584821 0.001279649 0.01208947      -6.370105 1.291674e-19
#AC092580.4 0.03730859 0.05112053 0.006003441 0.11743700      -3.090041 5.048838e-07
                   padj
#AAK1       1.294690e-28
#ABHD14A    1.708446e-07
#ABHD14B    4.636290e-05
#ABLIM1     1.540087e-02
#AC013264.2 1.950557e-15
#AC092580.4 7.254675e-03

# more examples 

# Comparing a condtion/condtions with different condtion/condtions (e.g. WT vs KO)
diff.res <- run.diff.exp(my.obj, de.by = "conditions", cond.1 = c("WT"), cond.2 = c("KO"))

# Comparing a cluster/clusters with different cluster/clusters (e.g. cluster 1 and 2 vs. 4)
diff.res <- run.diff.exp(my.obj, de.by = "clusters", cond.1 = c(1,4), cond.2 = c(2))

# Comparing a condtion/condtions with different condtion/condtions only in one/more cluster/clusters (e.g. cluster 1 WT vs cluster 1 KO)
diff.res <- run.diff.exp(my.obj, de.by = "clustBase.condComp", cond.1 = c("WT"), cond.2 = c("KO"), base.cond = 1)

# Comparing a cluster/clusters with different cluster/clusters only in one/more condition/conditions (e.g. cluster 1 vs cluster 2 but only the WT sample)
diff.res <- run.diff.exp(my.obj, de.by = "condBase.clustComp", cond.1 = c(1), cond.2 = c(2), base.cond = "WT")
```

- Volcano and MA plots 

```r
# Volcano Plot 
volcano.ma.plot(diff.res,
	sig.value = "pval",
	sig.line = 0.05,
	plot.type = "volcano",
	interactive = F)

# MA Plot
volcano.ma.plot(diff.res,
	sig.value = "pval",
	sig.line = 0.05,
	plot.type = "ma",
	interactive = F)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/volc_plot.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/MA_plot.png" width="400"/>      
</p>

 - Merging, resetting, renaming and removing clusters 
 
 ```r
# let's say you  want to merge cluster 3 and 2.
my.obj <- change.clust(my.obj, change.clust = 3, to.clust = 2)

# to reset to the original clusters run this.
my.obj <- change.clust(my.obj, clust.reset = T)

# you can also re-name the cluster numbers to cell types. Remember to reset after this so you can ran other analysis. 
my.obj <- change.clust(my.obj, change.clust = 7, to.clust = "B Cell")

# Let's say for what ever reason you want to remove acluster, to do so run this.
my.obj <- clust.rm(my.obj, clust.to.rm = 1)

# Remember that this would perminantly remove the data from all the slots in the object except frrom raw.data slot in the object. If you want to reset you need to start from the filtering cells step in the biginging of the analysis (using cell.filter function). 

# To re-position the cells run tSNE again 
my.obj <- run.tsne(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")

# Use this for plotting as you make the changes
cluster.plot(my.obj,
	cell.size = 1,
	plot.type = "tsne",
	cell.color = "black",
	back.col = "white",
	col.by = "clusters",
	cell.transparency = 0.5,
	clust.dim = 2,
	interactive = F)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_2D_a.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_2D_b.png" width="400"/>    
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_2D_c.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_2D_d.png" width="400"/>  
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_2D_e.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_2D_f.png" width="400"/>  
</p>

 - Cell gating 
 
  ```r
my.plot <- gene.plot(my.obj, gene = "GNLY", 
	plot.type = "scatterplot",
	clust.dim = 2,
	interactive = F)

cell.gating(my.obj, my.plot = my.plot)	

# or 

#my.plot <- cluster.plot(my.obj,
#	cell.size = 1,
#	cell.transparency = 0.5,
#	clust.dim = 2,
#	interactive = F)
 ```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/gate3.png" />
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/gate2.png" />    
</p>


After downloading the cell ids, use the following cammand to rename their cluster. 

```r
my.obj <- gate.to.clust(my.obj, my.gate = "cellGating.txt", to.clust = 10)
 ```

 - Pseudotime analysis
 
 ```r
MyGenes <- top.markers(marker.genes, topde = 50, min.base.mean = 0.2)
MyGenes <- unique(MyGenes)

pseudotime.tree(my.obj,
	marker.genes = MyGenes,
	type = "unrooted",
	clust.method = "complete")

# or 

pseudotime.tree(my.obj,
	marker.genes = MyGenes,
	type = "classic",
	clust.method = "complete")
	
pseudotime.tree(my.obj,
	marker.genes = MyGenes,
	type = "jitter",
	clust.method = "complete")	

 ```
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/bloodCells.jpg" />
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/pseudotime.tree_unrooted2.png" width="400" />
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/pseudotime.tree_unrooted.png" width="400" />
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/pseudotime.tree_classic.png" width="400" />
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/pseudotime.tree_jitter.png" width="400"/>
</p>


 - Pseudotime analysis using monocle
 
  ```r
library(monocle)

MyMTX <- my.obj@main.data
GeneAnno <- as.data.frame(row.names(MyMTX))
colnames(GeneAnno) <- "gene_short_name"
row.names(GeneAnno) <- GeneAnno$gene_short_name
cell.cluster <- (my.obj@best.clust)
Ha <- data.frame(do.call('rbind', strsplit(as.character(row.names(cell.cluster)),'_',fixed=TRUE)))[1]
clusts <- paste("cl.",as.character(cell.cluster$clusters),sep="")
cell.cluster <- cbind(cell.cluster,Ha,clusts)
colnames(cell.cluster) <- c("Clusts","iCellR.Conds","iCellR.Clusts")
Samp <- new("AnnotatedDataFrame", data = cell.cluster)
Anno <- new("AnnotatedDataFrame", data = GeneAnno)
my.monoc.obj <- newCellDataSet(as.matrix(MyMTX),phenoData = Samp, featureData = Anno)

## find disperesedgenes 
my.monoc.obj <- estimateSizeFactors(my.monoc.obj)
my.monoc.obj <- estimateDispersions(my.monoc.obj)
disp_table <- dispersionTable(my.monoc.obj)

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
my.monoc.obj <- setOrderingFilter(my.monoc.obj, unsup_clustering_genes$gene_id)

# tSNE
my.monoc.obj <- reduceDimension(my.monoc.obj, max_components = 2, num_dim = 10,reduction_method = 'tSNE', verbose = T)
# cluster 
my.monoc.obj <- clusterCells(my.monoc.obj, num_clusters = 10)

## plot conditions and clusters based on iCellR analysis 
A <- plot_cell_clusters(my.monoc.obj, 1, 2, color = "iCellR.Conds")
B <- plot_cell_clusters(my.monoc.obj, 1, 2, color = "iCellR.Clusts")

## plot clusters based monocle analysis 
C <- plot_cell_clusters(my.monoc.obj, 1, 2, color = "Cluster")

# get marker genes from iCellR analysis
MyGenes <- top.markers(marker.genes, topde = 30, min.base.mean = 0.2)
my.monoc.obj <- setOrderingFilter(my.monoc.obj, MyGenes)

my.monoc.obj <- reduceDimension(my.monoc.obj, max_components = 2,method = 'DDRTree')
# order cells 
my.monoc.obj <- orderCells(my.monoc.obj)

# plot based on iCellR analysis and marker genes from iCellR
D <- plot_cell_trajectory(my.monoc.obj, color_by = "iCellR.Clusts")

## heatmap genes from iCellR

plot_pseudotime_heatmap(my.monoc.obj[MyGenes,],
	cores = 1,
	cluster_rows = F,
	use_gene_short_name = T,
	show_rownames = T)
 ```
 
 <p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/13_monocol.png" />
	  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/14_monocol.png" />
</p>

# How to perform canonical correlation analysis (CCA) for sample alignment 

```r
# this function runs Seurat in the background but only works best if you have 2 samples. 
# If you have more than 2 samples it's best to run the commands explained in the next code block MNN.
# After CCA just replace the pca.data slot of iCellR with CCA results from Seurat. 

# We recommand running MNN alignment instead of CCA

library(iCellR)
library(Seurat)

# download sample 1
sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/sample1_for_CCA.tsv.gz"

download.file(url = sample.file.url, 
     destfile = "sample1_for_CCA.tsv.gz", 
     method = "auto")  


# download sample 2
sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/sample2_for_CCA.tsv.gz"

download.file(url = sample.file.url, 
     destfile = "sample2_for_CCA.tsv.gz", 
     method = "auto")  

# Read both samples 
S1 <- read.table("sample1_for_CCA.tsv.gz")
head(S1)[1:5]

S2 <- read.table("sample2_for_CCA.tsv.gz")
head(S2)[1:5]

# aggregate both samples  
my.data <- data.aggregation(samples = c("S1","S2"), condition.names = c("S1","S2"))

# make object
my.obj <- make.obj(my.data)

# QC
my.obj <- qc.stats(my.obj,
s.phase.genes = s.phase, 
g2m.phase.genes = g2m.phase)

# filter
my.obj <- cell.filter(my.obj)

## CCA 

require(devtools)
install_version("Seurat", version = "2.3.4", repos = "http://cran.us.r-project.org")
library(Seurat)

# this function runs Seurat in the background but only works best if you have 2 samples. 
# If you have more than 2 samples it's best to run the commands explained in the next code block MNN.
# After CCA just replace the pca.data slot of iCellR with CCA results from Seurat. 

my.obj <- run.cca(my.obj,
	top.vari.genes = 1000,
	cc.number = 30,
	dims.align = 1:20,
	normalize.data = TRUE,
	scale.data = TRUE,
	normalization.method = "LogNormalize",
	scale.factor = 10000,
	display.progress = TRUE)

################ Normalize, gene stats and gene model for PCA
my.obj <- norm.data(my.obj, norm.method = "ranked.glsf", top.rank = 500) 
my.obj <- gene.stats(my.obj, which.data = "main.data")
my.obj <- make.gene.model(my.obj)

###### See data without CCA
my.obj <- run.pca(my.obj, method = "gene.model", gene.list = my.obj@gene.model,data.type = "main",batch.norm = F)

my.obj <- run.umap(my.obj, dims = 1:10, method = "umap-learn") 

UMAP_NoCCA <- cluster.plot(my.obj,plot.type = "umap",cell.color = "black",col.by = "conditions",cell.transparency = 0.5,interactive = F)

PCA <- cluster.plot(my.obj,plot.type = "pca",cell.color = "black",col.by = "conditions",cell.transparency = 0.5,interactive = F)

##### See data with CCA
# replace PCA with CCA
attributes(my.obj)$pca.data <- my.obj@cca.data

my.obj <- run.umap(my.obj, dims = 1:10, method = "umap-learn") 

UMAP_CCA <- cluster.plot(my.obj,plot.type = "umap",cell.color = "black",col.by = "conditions",cell.transparency = 0.5,interactive = F)

Aligned_CCA <- cluster.plot(my.obj,plot.type = "pca",cell.color = "black",col.by = "conditions",cell.transparency = 0.5,interactive = F)

# plot
library(gridExtra)
png('Compare.png', width = 8, height = 8, units = 'in', res = 300)
grid.arrange(UMAP_NoCCA,PCA,UMAP_CCA,Aligned_CCA)
dev.off()

```
Before and After CCA

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/Compare.png" />
</p>


# How to perform mutual nearest neighbor (MNN) sample alignment 

```r
library(iCellR)

# download sample 1
sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/sample1_for_CCA.tsv.gz"

download.file(url = sample.file.url, 
     destfile = "sample1_for_CCA.tsv.gz", 
     method = "auto")  


# download sample 2
sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/sample2_for_CCA.tsv.gz"

download.file(url = sample.file.url, 
     destfile = "sample2_for_CCA.tsv.gz", 
     method = "auto")  

# Read both samples 
S1 <- read.table("sample1_for_CCA.tsv.gz")
head(S1)[1:5]

S2 <- read.table("sample2_for_CCA.tsv.gz")
head(S2)[1:5]

# aggregate both samples  
my.data <- data.aggregation(samples = c("S1","S2"), condition.names = c("S1","S2"))

# make object
my.obj <- make.obj(my.data)

# QC
my.obj <- qc.stats(my.obj,
s.phase.genes = s.phase, 
g2m.phase.genes = g2m.phase)

# filter
my.obj <- cell.filter(my.obj)
my.obj <- gene.stats(my.obj, which.data = "main.data")

my.obj <- make.gene.model(my.obj, my.out.put = "data",
	dispersion.limit = 1.5,
	base.mean.rank = 500,
	no.mito.model = T,
	mark.mito = T,
	interactive = F,
	no.cell.cycle = T,
	out.name = "gene.model")
	
###### Run MNN 
# This would automatically run all the samples in your experiment 

library(scran)
my.obj <- run.mnn(my.obj,
    top.rank = 500,
    k=20,
    d=50)

# detach the scran pacakge after MNN as it masks some of the functions 
detach("package:scran", unload=TRUE)

# or 
#my.obj <- run.mnn(my.obj,
#    method = "gene.model",
#    gene.list = my.obj@gene.model,
#    k=20,
#    d=50)

# By running MNN alignment you replace the PCA data slot and there is no need to run PCA. 
# If you run PCA, MNN results will be replaced. (Do only if you want to see the results before and after MNN)

# detach scran
detach("package:scran", unload=TRUE)

# normaliza the main data for iCellR analyses
my.obj <- norm.data(my.obj, norm.method = "ranked.glsf", top.rank = 500)

## run tSNE 
my.obj <- run.pc.tsne(my.obj, dims = 1:10)
my.obj <- run.umap(my.obj, dims = 1:10, method = "umap-learn")

cluster.plot(my.obj,plot.type = "tsne",col.by = "conditions",interactive = F)

cluster.plot(my.obj,plot.type = "umap",col.by = "conditions",interactive = F)
```

Before and After MNN analysis

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/MNN.png" />
</p>


# How to perform Seurat's anchor (integration) sample alignment 

```r
library(iCellR)

# download sample 1
sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/sample1_for_CCA.tsv.gz"

download.file(url = sample.file.url, 
     destfile = "sample1_for_CCA.tsv.gz", 
     method = "auto")  


# download sample 2
sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/sample2_for_CCA.tsv.gz"

download.file(url = sample.file.url, 
     destfile = "sample2_for_CCA.tsv.gz", 
     method = "auto")  

# Read both samples 
S1 <- read.table("sample1_for_CCA.tsv.gz")
head(S1)[1:5]

S2 <- read.table("sample2_for_CCA.tsv.gz")
head(S2)[1:5]

# aggregate both samples  
my.data <- data.aggregation(samples = c("S1","S2"), condition.names = c("S1","S2"))

# make object
my.obj <- make.obj(my.data)

# QC
my.obj <- qc.stats(my.obj,
s.phase.genes = s.phase, 
g2m.phase.genes = g2m.phase)

# filter
my.obj <- cell.filter(my.obj)
my.obj <- gene.stats(my.obj, which.data = "main.data")

my.obj <- make.gene.model(my.obj, my.out.put = "data",
	dispersion.limit = 1.5,
	base.mean.rank = 500,
	no.mito.model = T,
	mark.mito = T,
	interactive = F,
	no.cell.cycle = T,
	out.name = "gene.model")
	
###### Run Anchor 
# This would automatically run all the samples in your experiment 

library(Seurat)
my.obj <- run.anchor(my.obj,
    top.rank = 500,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    selection.method = "vst",
    nfeatures = 2000,
    dims = 1:20)

# By running Anchor alignment you replace the PCA data slot and there is no need to run PCA. 
# If you run PCA, Anchor results will be replaced. (Do only if you want to see the results before and after MNN)

# normaliza the main data for iCellR analyses
my.obj <- norm.data(my.obj, norm.method = "ranked.glsf", top.rank = 500)

## run tSNE 
my.obj <- run.pc.tsne(my.obj, dims = 1:10)
my.obj <- run.umap(my.obj, dims = 1:10, method = "umap-learn")

cluster.plot(my.obj,plot.type = "tsne",col.by = "conditions",interactive = F)

cluster.plot(my.obj,plot.type = "umap",col.by = "conditions",interactive = F)
```

Before and After Anchor analysis

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/Anchor.png" />
</p>

# How to demultiplex with hashtag oligos (HTOs)

```r
# Read an example file
 
my.hto <- read.table(file = system.file('extdata', 'dense_umis.tsv', package = 'iCellR'), as.is = TRUE)
 
# see the head of the file for the first few columns
head(my.hto)[1:3]
#                         TGACAACAGGGCTCTC AAGGAGCGTCATTAGC AGTGAGGAGACTGTAA
#Hashtag1-GTCAACTCTTTAGCG                3                7                7
#Hashtag2-TGATGGCCTATTGGG               18               24             1761
#Hashtag3-TTCCGCCTCTCTTTG                7                8                5
#Hashtag4-AGTAAGTTCAGCGTA                0                0                0
#Hashtag5-AAGTATCGTTTCGCA              890                2               11
#Hashtag7-TGTCTTTCCTGCCAG                5                3                3
 
# run annotation
htos <- hto.anno(hto.data = my.hto, cov.thr = 10, assignment.thr = 80)

head(htos)
#                 Hashtag1-GTCAACTCTTTAGCG Hashtag2-TGATGGCCTATTGGG
#TGACAACAGGGCTCTC                        3                       18
#AAGGAGCGTCATTAGC                        7                       24
#AGTGAGGAGACTGTAA                        7                     1761
#ATCCACCCATGTTCCC                      753                       20
#AAACGGGCAGGACCCT                      728                       24
#ATGTGTGAGTCTTGCA                        4                       25
#                 Hashtag3-TTCCGCCTCTCTTTG Hashtag4-AGTAAGTTCAGCGTA
#TGACAACAGGGCTCTC                        7                        0
#AAGGAGCGTCATTAGC                        8                        0
#AGTGAGGAGACTGTAA                        5                        0
#ATCCACCCATGTTCCC                        3                        0
#AAACGGGCAGGACCCT                        3                        0
#ATGTGTGAGTCTTGCA                      370                        0
#                 Hashtag5-AAGTATCGTTTCGCA Hashtag7-TGTCTTTCCTGCCAG unmapped
#TGACAACAGGGCTCTC                      890                        5       17
#AAGGAGCGTCATTAGC                        2                        3        3
#AGTGAGGAGACTGTAA                       11                        3       87
#ATCCACCCATGTTCCC                        5                        6       18
#AAACGGGCAGGACCCT                        9                        3       16
#ATGTGTGAGTCTTGCA                        9                     1011       25
#                    assignment.annotation percent.match coverage low.cov
#TGACAACAGGGCTCTC Hashtag5-AAGTATCGTTTCGCA      94.68085      940   FALSE
#AAGGAGCGTCATTAGC Hashtag2-TGATGGCCTATTGGG      51.06383       47    TRUE
#AGTGAGGAGACTGTAA Hashtag2-TGATGGCCTATTGGG      93.97012     1874   FALSE
#ATCCACCCATGTTCCC Hashtag1-GTCAACTCTTTAGCG      93.54037      805   FALSE
#AAACGGGCAGGACCCT Hashtag1-GTCAACTCTTTAGCG      92.97573      783   FALSE
#ATGTGTGAGTCTTGCA Hashtag7-TGTCTTTCCTGCCAG      70.01385     1444   FALSE
#                 assignment.threshold
#TGACAACAGGGCTCTC      good.assignment
#AAGGAGCGTCATTAGC               unsure
#AGTGAGGAGACTGTAA      good.assignment
#ATCCACCCATGTTCCC      good.assignment
#AAACGGGCAGGACCCT      good.assignment
#ATGTGTGAGTCTTGCA               unsure

# plot

A = ggplot(htos, aes(assignment.annotation,percent.match)) +
	geom_jitter(alpha = 0.25, color = "blue") +
	geom_boxplot(alpha = 0.5) + 
	theme_bw() + 
	theme(axis.text.x=element_text(angle=90))

B = ggplot(htos, aes(low.cov,percent.match)) +
	geom_jitter(alpha = 0.25, color = "blue") +
	geom_boxplot(alpha = 0.5) + 
	theme_bw() + 
	theme(axis.text.x=element_text(angle=90))

library(gridExtra)
png('HTOs.png', width = 8, height = 8, units = 'in', res = 300)
grid.arrange(A,B,ncol=2)
dev.off()
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/data/myHTOS.png" />
</p>

 - Filtering HTOs and merging the samples
 
 ```r
 # let's say you decided filtering based on 80%
 dim(htos)
 # [1] 1500   12
 htos <- subset(htos,htos$percent.match > 80)
 dim(htos)
 # [1] 1073   12 
 
 # Take the cell IDs from Hashtag1
 sample1 <- row.names(subset(htos,htos$assignment.annotation == "Hashtag1-GTCAACTCTTTAGCG"))
 
 head(sample1)
# [1] "ATCCACCCATGTTCCC" "AAACGGGCAGGACCCT" "TTCTACATCCTCATTA" "GGTATTGTCCTATGTT"
# [5] "GTCGTAATCTTACCTA" "ACAGCCGGTTGGGACA"

length(sample1)
# [1] 213
# in this case you have 213 cells in sample 1 (Hashtag1)

# Take the cell IDs from Hashtag2
sample2 <- row.names(subset(htos,htos$assignment.annotation == "Hashtag2-TGATGGCCTATTGGG"))

# now read your RNA data 
# example:
RNA.data <- load10x("YOUR/data/filtered_gene_bc_matrices/hg19/")

head(RNA.data)[1:2]
#         AAACATAAAACCAG CCCCATTGAGCTAA
#A1BG.AS1                   0                   0
#BCLA                       0                   0
#A2M                        0                   0
#GATA1                      0                   0

# NOTE: the RNA data has the cell IDs in the same format as HTOs 
# "AAACATAAAACCAG" "CCCCATTGAGCTAA" ... 
# Not "AAACATAAAACCAG.1" "CCCCATTGAGCTAA.1" ... 

# demultiplex RNA data 
# Take RNA-Seq data with the cell IDs from Hashtag1
sample1.rna <- RNA.data[ , which(names(RNA.data) %in% sample1)]

# Take RNA-Seq data with the cell IDs from Hashtag2
sample2.rna <- RNA.data[ , which(names(RNA.data) %in% sample2)]

# aggregate (merge the 2 or more samples after demultiplexing)

my.data <- data.aggregation(samples = c("sample1.rna","sample2.rna"), 
	condition.names = c("S1","S2"))
	
# make iCellR object	
my.obj <- make.obj(my.data)

# The rest is as above :)
 ```


# How to analyze CITE-seq data using iCellR

 - Download test samples
 
 ```r
 sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/CITE-Seq_sample_RNA.tsv.gz"

# download RNA file

download.file(url = sample.file.url, 
     destfile = "CITE-Seq_sample_RNA.tsv.gz", 
     method = "auto")  

sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/CITE-Seq_sample_ADT.tsv.gz"

# download ADT file

download.file(url = sample.file.url, 
     destfile = "CITE-Seq_sample_ADT.tsv.gz", 
     method = "auto")  
 ```
 
  - Read the files and make your object
  
  ```r
  # Read RNA file
 rna.data <- read.delim("CITE-Seq_sample_RNA.tsv.gz",header=TRUE)
 
  # see the head 
 head(rna.data)[1:3]
#          CTGTTTACACCGCTAG CTCTACGGTGTGGCTC AGCAGCCAGGCTCATT
#A1BG                    0                0                0
#A1BG-AS1                0                0                0
#A1CF                    0                0                0
#A2M                     0                0                0
#A2M-AS1                 0                0                0
#A2ML1                   0                0                0
  
 # Read ADT file
 adt.data <- read.delim("CITE-Seq_sample_ADT.tsv.gz",header=TRUE)
 
 # see the head 
 head(adt.data)[1:3]
#        CTGTTTACACCGCTAG CTCTACGGTGTGGCTC AGCAGCCAGGCTCATT
#CD3                  60               52               89
#CD4                  72               49              112
#CD8                  76               59               61
#CD45RA              575             3943              682
#CD56                 64               68               87
#CD16                161              107              117
 
# if you had multiple sample use the data.aggregation function for both RNA and ADT data. 

# make iCellR object
my.obj <- make.obj(rna.data)

# check object
my.obj
###################################
,--. ,-----.       ,--.,--.,------.
`--''  .--./ ,---. |  ||  ||  .--. '
,--.|  |    | .-. :|  ||  ||  '--'.'
|  |'  '--'\   --. |  ||  ||  |
`--' `-----' `----'`--'`--'`--' '--'
###################################
An object of class iCellR version: 1.1.4
Raw/original data dimentions (rows,columns): 20501,8617
Data conditions: no conditions/single sample
Row names: A1BG,A1BG-AS1,A1CF ...
Columns names: CTGTTTACACCGCTAG,CTCTACGGTGTGGCTC,AGCAGCCAGGCTCATT ...
###################################
   QC stats performed:FALSE, PCA performed:FALSE, CCA performed:FALSE
   Clustering performed:FALSE, Number of clusters:0
   tSNE performed:FALSE, UMAP performed:FALSE, DiffMap performed:FALSE
   Main data dimentions (rows,columns):0,0
   Normalization factors:,...
   Imputed data dimentions (rows,columns):0,0
############## scVDJ-Seq ###########
VDJ data dimentions (rows,columns):0,0
############## CITE-Seq ############
   ADT raw data dimentions (rows,columns):0,0
   ADT main data dimentions (rows,columns):0,0
   ADT columns names:...
   ADT row names:...
########### iCellR object ##########
```

- add ADT data

```r
my.obj <- add.adt(my.obj, adt.data = adt.data)

# check too see
 my.obj
###################################
,--. ,-----.       ,--.,--.,------.
`--''  .--./ ,---. |  ||  ||  .--. '
,--.|  |    | .-. :|  ||  ||  '--'.'
|  |'  '--'\   --. |  ||  ||  |
`--' `-----' `----'`--'`--'`--' '--'
###################################
An object of class iCellR version: 1.1.4
Raw/original data dimentions (rows,columns): 20501,8617
Data conditions: no conditions/single sample
Row names: A1BG,A1BG-AS1,A1CF ...
Columns names: CTGTTTACACCGCTAG,CTCTACGGTGTGGCTC,AGCAGCCAGGCTCATT ...
###################################
   QC stats performed:FALSE, PCA performed:FALSE, CCA performed:FALSE
   Clustering performed:FALSE, Number of clusters:0
   tSNE performed:FALSE, UMAP performed:FALSE, DiffMap performed:FALSE
   Main data dimentions (rows,columns):0,0
   Normalization factors:,...
   Imputed data dimentions (rows,columns):0,0
############## scVDJ-Seq ###########
VDJ data dimentions (rows,columns):0,0
############## CITE-Seq ############
-   ADT raw data dimentions (rows,columns):10,8617
   ADT main data dimentions (rows,columns):0,0
   ADT columns names:...
   ADT row names:...
########### iCellR object ##########
  ```
- QC, filter, normalize, merge ADT and RNA data, run PCA and UMAP

```r
# QC
my.obj <- qc.stats(my.obj,
	s.phase.genes = s.phase, 
	g2m.phase.genes = g2m.phase)

# plot as mentioned above

# filter 
my.obj <- cell.filter(my.obj,
	min.mito = 0,
	max.mito = 0.07 ,
	min.genes = 500,
	max.genes = 4000,
	min.umis = 0,
	max.umis = Inf)

# normalize RNA
my.obj <- norm.data(my.obj, norm.method = "ranked.glsf", top.rank = 500) 

# normalize ADT
my.obj <- norm.adt(my.obj)

# gene stats
my.obj <- gene.stats(my.obj, which.data = "main.data")

# find genes for PCA
my.obj <- make.gene.model(my.obj, my.out.put = "data",
	dispersion.limit = 1.5, 
	base.mean.rank = 500, 
	no.mito.model = T, 
	mark.mito = T, 
	interactive = F,
	no.cell.cycle = T,
	out.name = "gene.model")

# merge RNA and ADT data
my.obj <- adt.rna.merge(my.obj, adt.data = "main")

# run PCA and the rest is as above

my.obj <- run.pca(my.obj, method = "gene.model", gene.list = my.obj@gene.model,data.type = "main",batch.norm = F)

# 2 pass PCA 
my.obj <- find.dim.genes(my.obj, dims = 1:20,top.pos = 20, top.neg = 20)
# second round PC
my.obj <- run.pca(my.obj, method = "gene.model", gene.list = my.obj@gene.model,data.type = "main",batch.norm = F)

my.obj <- run.umap(my.obj, dims = 1:10, method = "umap-learn")

# check your object 
my.obj
###################################
,--. ,-----.       ,--.,--.,------.
`--''  .--./ ,---. |  ||  ||  .--. '
,--.|  |    | .-. :|  ||  ||  '--'.'
|  |'  '--'\   --. |  ||  ||  |
`--' `-----' `----'`--'`--'`--' '--'
###################################
An object of class iCellR version: 1.1.4
Raw/original data dimentions (rows,columns): 20501,8617
Data conditions: no conditions/single sample
Row names: A1BG,A1BG-AS1,A1CF ...
Columns names: CTGTTTACACCGCTAG,CTCTACGGTGTGGCTC,AGCAGCCAGGCTCATT ...
###################################
   QC stats performed:TRUE, PCA performed:TRUE, CCA performed:FALSE
   Clustering performed:TRUE, Number of clusters:14
   tSNE performed:FALSE, UMAP performed:TRUE, DiffMap performed:FALSE
   Main data dimentions (rows,columns):20511,8305
   Normalization factors:8.448547776071,...
   Imputed data dimentions (rows,columns):0,0
############## scVDJ-Seq ###########
VDJ data dimentions (rows,columns):0,0
############## CITE-Seq ############
   ADT raw data dimentions (rows,columns):10,8617
   ADT main data dimentions (rows,columns):10,8617
   ADT columns names:CTGTTTACACCGCTAG...
   ADT row names:ADT_CD3...
########### iCellR object ##########
```

- plot 

```r
# find ADT gene names 
grep("^ADT_", rownames(my.obj@main.data),value=T)
# [1] "ADT_CD3"    "ADT_CD4"    "ADT_CD8"    "ADT_CD45RA" "ADT_CD56"
# [6] "ADT_CD16"   "ADT_CD11c"  "ADT_CD14"   "ADT_CD19"   "ADT_CD34"

A = gene.plot(my.obj, 
	gene = "ADT_CD3",
	plot.data.type = "umap",
	interactive = F,
	cell.transparency = 0.5)

B = gene.plot(my.obj, 
	gene = "CD3E",
	plot.data.type = "umap",
	interactive = F,
	cell.transparency = 0.5)

C = gene.plot(my.obj, 
	gene = "ADT_CD16",
	plot.data.type = "umap",
	interactive = F,
	cell.transparency = 0.5)

D = gene.plot(my.obj, 
	gene = "FCGR3A",
	plot.data.type = "umap",
	interactive = F,
	cell.transparency = 0.5)
		
library(gridExtra)
grid.arrange(A,B,C,D)
```


<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/NotScaled.png" />
</p>

# How to analyze scVDJ-seq data using iCellR

Here is an example of how to add VDJ data.

 ```r
# first prepare the files. 
# this function would filter the files, calculate clonotype frequencies and proportions and add conditions to the cell ids.
my.vdj.1 <- prep.vdj(vdj.data = "all_contig_annotations.csv", cond.name = "WT")
my.vdj.2 <- prep.vdj(vdj.data = "all_contig_annotations.csv", cond.name = "KO")
my.vdj.3 <- prep.vdj(vdj.data = "all_contig_annotations.csv", cond.name = "Ctrl")

# concatenate all the conditions
my.vdj.data <- rbind(my.vdj.1, my.vdj.2, my.vdj.3)

# see head of the file
head(my.vdj.data)
#  raw_clonotype_id               barcode is_cell                   contig_id
#1       clonotype1 WT_AAACCTGAGCTAACTC-1    True AAACCTGAGCTAACTC-1_contig_1
#2       clonotype1 WT_AAACCTGAGCTAACTC-1    True AAACCTGAGCTAACTC-1_contig_2
#3       clonotype1 WT_AGTTGGTTCTCGCATC-1    True AGTTGGTTCTCGCATC-1_contig_3
#4       clonotype1 WT_TGACAACCAACTGCTA-1    True TGACAACCAACTGCTA-1_contig_1
#5       clonotype1 WT_TGTCCCAGTCAAACTC-1    True TGTCCCAGTCAAACTC-1_contig_1
#6       clonotype1 WT_TGTCCCAGTCAAACTC-1    True TGTCCCAGTCAAACTC-1_contig_2
#  high_confidence length chain  v_gene d_gene  j_gene c_gene full_length
#1            True    693   TRA TRAV8-1   None  TRAJ21   TRAC        True
#2            True    744   TRB  TRBV28  TRBD1 TRBJ2-1  TRBC2        True
#3            True    647   TRA TRAV8-1   None  TRAJ21   TRAC        True
#4            True    508   TRB  TRBV28  TRBD1 TRBJ2-1  TRBC2        True
#5            True    660   TRA TRAV8-1   None  TRAJ21   TRAC        True
#6            True    770   TRB  TRBV28  TRBD1 TRBJ2-1  TRBC2        True
#  productive             cdr3                                          cdr3_nt
#1       True      CAVKDFNKFYF                TGTGCCGTGAAAGACTTCAACAAATTTTACTTT
#2       True CASSLFSGTGTNEQFF TGTGCCAGCAGTTTATTTTCCGGGACAGGGACGAATGAGCAGTTCTTC
#3       True      CAVKDFNKFYF                TGTGCCGTGAAAGACTTCAACAAATTTTACTTT
#4       True CASSLFSGTGTNEQFF TGTGCCAGCAGTTTATTTTCCGGGACAGGGACGAATGAGCAGTTCTTC
#5       True      CAVKDFNKFYF                TGTGCCGTGAAAGACTTCAACAAATTTTACTTT
#6       True CASSLFSGTGTNEQFF TGTGCCAGCAGTTTATTTTCCGGGACAGGGACGAATGAGCAGTTCTTC
#  reads umis       raw_consensus_id my.raw_clonotype_id clonotype.Freq
#1  1241    2 clonotype1_consensus_1          clonotype1            120
#2  2400    4 clonotype1_consensus_2          clonotype1            120
#3  1090    2 clonotype1_consensus_1          clonotype1            120
#4  2455    4 clonotype1_consensus_2          clonotype1            120
#5  1346    2 clonotype1_consensus_1          clonotype1            120
#6  3073    8 clonotype1_consensus_2          clonotype1            120
#  proportion total.colonotype
#1 0.04098361             1292
#2 0.04098361             1292
#3 0.04098361             1292
#4 0.04098361             1292
#5 0.04098361             1292
#6 0.04098361             1292

# add it to iCellR object
add.vdj(my.obj, vdj.data = my.vdj.data)
 ```

# How to analyze large bulk RNA-Seq data (TCGA)

In this example the samples are normalized using DESeq2 so no noramalizaion is needed.

```r
sample.file.url = "https://genome.med.nyu.edu/results/external/iCellR/data/TCGA_sample_Normalized_data.tsv.gz"

download.file(url = sample.file.url, 
     destfile = "TCGA_sample_Normalized_data.tsv.gz", 
     method = "auto")  

TCGA.data <- read.table("TCGA_sample_Normalized_data.tsv.gz")
head(TCGA.data)[1:3]
#         Basal_TCGA.A1.A0SK.txt Basal_TCGA.A1.A0SP.txt Basal_TCGA.A2.A04P.txt
#TSPAN6                5823.4300            4318.034382            5265.733258
#TNMD                     0.0000               6.049079               6.763079
#DPM1                  3248.1536            2528.515113            1183.538813
#SCYL3                 1059.7135             965.836315            1109.144945
#C1orf112              1251.3155            1070.687022             485.589067
#FGR                    106.2438             933.574559             512.641383

library(iCellR)
my.obj <- make.obj(TCGA.data)

my.obj@main.data <- my.obj@raw.data

my.obj
###################################
,--. ,-----.       ,--.,--.,------.
`--''  .--./ ,---. |  ||  ||  .--. '
,--.|  |    | .-. :|  ||  ||  '--'.'
|  |'  '--'\   --. |  ||  ||  |
`--' `-----' `----'`--'`--'`--' '--'
###################################
An object of class iCellR version: 1.2.4
Raw/original data dimentions (rows,columns): 69797,882
Data conditions in raw data: Basal,Her2,LumA,LumB,Normal (131,64,404,170,113)
Row names: TSPAN6,TNMD,DPM1 ...
Columns names: Basal_TCGA.A1.A0SK.txt,Basal_TCGA.A1.A0SP.txt,Basal_TCGA.A2.A04P.txt ...
###################################
   QC stats performed:FALSE, PCA performed:FALSE, CCA performed:FALSE
   Clustering performed:FALSE, Number of clusters:0
   tSNE performed:FALSE, UMAP performed:FALSE, DiffMap performed:FALSE
   Main data dimentions (rows,columns):69797,882
   Normalization factors:,...
   Imputed data dimentions (rows,columns):0,0
############## scVDJ-Seq ###########
VDJ data dimentions (rows,columns):0,0
############## CITE-Seq ############
   ADT raw data dimentions (rows,columns):0,0
   ADT main data dimentions (rows,columns):0,0
   ADT columns names:...
   ADT row names:...
########### iCellR object ##########


my.obj <- run.pca(my.obj)

my.obj <- run.clustering(my.obj, 
	clust.method = "kmeans", 
	dist.method = "euclidean",
	index.method = "silhouette",
	max.clust =25,
	min.clust = 2,
	dims = 1:10)

my.obj <- run.pc.tsne(my.obj, dims = 1:10)
my.obj <- run.umap(my.obj, dims = 1:10, method = "umap-learn") 

cluster.plot(my.obj,plot.type = "pca",cell.color = "black",col.by = "conditions",cell.transparency = 0.5,interactive = F)
cluster.plot(my.obj,plot.type = "umap",cell.color = "black",col.by = "conditions",cell.transparency = 0.5,interactive = F)
cluster.plot(my.obj,plot.type = "tsne",cell.color = "black",col.by = "conditions",cell.transparency = 0.5,interactive = F)
cluster.plot(my.obj,plot.type = "umap",cell.color = "black",cell.transparency = 1,interactive = F)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/TCGA.png" />
</p>


```r
sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux

Matrix products: default
BLAS: /gpfs/share/apps/R/3.5.1/lib64/R/lib/libRblas.so
LAPACK: /gpfs/share/apps/R/3.5.1/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] iCellR_1.2.2  plotly_4.9.0  ggplot2_3.2.1

loaded via a namespace (and not attached):
 [1] ggrepel_0.8.1        Rcpp_1.0.2           ape_5.3
 [4] lattice_0.20-38      tidyr_1.0.0          assertthat_0.2.1
 [7] zeallot_0.1.0        digest_0.6.22        mime_0.7
[10] R6_2.4.0             plyr_1.8.4           backports_1.1.5
[13] acepack_1.4.1        httr_1.4.1           pillar_1.4.2
[16] rlang_0.4.0          lazyeval_0.2.2       rstudioapi_0.10
[19] data.table_1.12.6    rpart_4.1-15         Matrix_1.2-17
[22] checkmate_1.9.4      reticulate_1.13      splines_3.5.1
[25] Rtsne_0.15           stringr_1.4.0        foreign_0.8-72
[28] htmlwidgets_1.5.1    pheatmap_1.0.12      munsell_0.5.0
[31] umap_0.2.3.1         shiny_1.4.0          compiler_3.5.1
[34] httpuv_1.5.2         xfun_0.10            pkgconfig_2.0.3
[37] askpass_1.1          base64enc_0.1-3      htmltools_0.4.0
[40] nnet_7.3-12          openssl_1.4.1        tidyselect_0.2.5
[43] htmlTable_1.13.2     tibble_2.1.3         gridExtra_2.3
[46] Hmisc_4.2-0          reshape_0.8.8        viridisLite_0.3.0
[49] ggpubr_0.2.3         crayon_1.3.4         dplyr_0.8.3
[52] withr_2.1.2          later_1.0.0          MASS_7.3-51.4
[55] grid_3.5.1           NbClust_3.0          nlme_3.1-141
[58] jsonlite_1.6         xtable_1.8-4         gtable_0.3.0
[61] lifecycle_0.1.0      magrittr_1.5         scales_1.0.0
[64] stringi_1.4.3        ggsignif_0.6.0       promises_1.1.0
[67] scatterplot3d_0.3-41 latticeExtra_0.6-28  ggdendro_0.1-20
[70] vctrs_0.2.0          Formula_1.2-3        RColorBrewer_1.1-2
[73] tools_3.5.1          glue_1.3.1           purrr_0.3.3
[76] parallel_3.5.1       fastmap_1.0.1        survival_2.44-1.1
[79] colorspace_1.4-1     cluster_2.1.0        knitr_1.25
```
