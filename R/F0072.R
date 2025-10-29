#' Create bubble heatmaps for genes in clusters or conditions.
#'
#' This function takes an object of class iCellR and genes and provides a heatmap.
#' @param x A data frame containing gene counts for cells.
#' @param gene A set of gene names to be heatmapped.
#' @param conds.to.plot Choose the conditions you want to see in the plot, default = NULL (all conditions).
#' @param heat.colors Colors for heatmap, default = c("blue" ,"white", "red").
#' @param interactive If TRUE an html interactive file will be made, default = TRUE.
#' @param out.name Output name for html file if interactive = TRUE, default = "plot".
#' @param data.type Choose from "main", "atac", atac.imputed and "imputed", default = "main".
#' @param min.scale Set a minimum color scale, default = -2.5.
#' @param max.scale Set a maximum color scale, default = 2.5.
#' @param colour Set color to "Percent.Expressed", or "Expression", , default = "Expression".
#' @param size Set size to "Percent.Expressed", or "Expression", , default = "Percent.Expressed".
#' @param write.data Write export the data used for the plot plot, default = TFALSE.
#' @return An object of class iCellR
#' @import pheatmap
#' @importFrom reshape melt
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly ggplotly layout plot_ly
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @importFrom methods new
#' @importFrom stats aggregate as.dendrogram cor cor.test dist hclust p.adjust prcomp quantile sd t.test
#' @importFrom utils capture.output packageVersion read.table write.table
#' @importFrom graphics legend par plot
#' @importFrom ggplot2 ggplot geom_segment geom_violin guide_colorbar guide_legend guides scale_color_discrete scale_colour_gradient scale_fill_gradient2 scale_x_continuous scale_y_continuous scale_y_discrete stat_summary coord_polar element_rect element_text element_blank facet_wrap scale_color_manual geom_hline geom_jitter geom_vline ylab xlab ggtitle theme_bw aes theme geom_bar geom_point geom_boxplot geom_errorbar position_dodge geom_tile geom_density geom_line
#' @export
bubble.gg.plot <- function (x = NULL,
                          gene = "NULL",
                          data.type = "main",
                          conds.to.plot = NULL,
                          min.scale = -2.5,
                          max.scale = 2.5,
                          interactive = TRUE,
                          write.data = FALSE,
                          colour = "Expression",
                          size = "Percent.Expressed",
                          out.name = "plot",
                          heat.colors = c("blue","red")) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ## get main data

##################################
DATA <- x@best.clust
MY.conds <- as.data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
MY.conds <- as.character(as.matrix(MY.conds))
DATA <- as.data.frame(cbind(DATA,MY.conds))
#######
if (!is.null(conds.to.plot)) {
  DATA <- subset(DATA, DATA$MY.conds %in% conds.to.plot)
  message(paste(" "))
  message(paste("#############"))
  message(paste(" Averaging gene expression for all clusters in:",paste(conds.to.plot,collapse = ","),"..."))
  message(paste("#############"))
  message(paste(" "))
}
###########
# get data
sampleCondition <- DATA$clusters
conditions <- sort(unique(sampleCondition))
DATA1 <- DATA

## get main data
if (data.type == "main") {
  Table <- x@main.data
}
if (data.type == "imputed") {
  Table <- x@imputed.data
}
if (data.type == "atac") {
  Table <- x@atac.main
}
if (data.type == "atac.imputed") {
  Table <- x@atac.imputed
}

#################################### get the genes
Table  <- Table [gene, colnames(Table),drop = FALSE]

#  Table = x@main.data
datalist <- list()
###
for(i in conditions){
  IDs <- rownames(subset(DATA1, sampleCondition == i))
  DATA <- Table[ , which(names(Table) %in% IDs)]
  DATA <- as.matrix(DATA)
  DATA <- apply(DATA, 1, function(DATA) {mean(DATA)})
  DATA <- as.data.frame(DATA)
  Name=paste("meanExp_cluster",i,".txt",sep="_")
  NameCol=paste("cluster",i,sep="_")
  colnames(DATA) <- NameCol
  DATA <- cbind(gene = rownames(DATA), DATA)
  rownames(DATA) <- NULL
  eval(call("<-", as.name(NameCol), DATA))
}
########### reduce
filenames <- ls(pattern="cluster_")
filenames <- filenames[order(nchar(filenames))]
datalist <- mget(filenames)
MeanExpForClusters <- Reduce(function(x,y) {merge(x,y)}, datalist)

################### percent expressed
for(i in conditions){
  IDs <- rownames(subset(DATA1, sampleCondition == i))
  DATA <- Table[ , which(names(Table) %in% IDs)]
  DATA <- as.matrix(DATA)
#
  DATA <- (rowSums(DATA > 0) / dim(DATA)[2])*100
#
  DATA <- as.data.frame(DATA)
  Name=paste("meanExp_cluster",i,".txt",sep="_")
  NameCol=paste("cluster",i,sep="_")
  colnames(DATA) <- NameCol
  DATA <- cbind(gene = rownames(DATA), DATA)
  rownames(DATA) <- NULL
  eval(call("<-", as.name(NameCol), DATA))
}
########### reduce
filenames <- ls(pattern="cluster_")
filenames <- filenames[order(nchar(filenames))]
datalist <- mget(filenames)
MeanExpForClusters1 <- Reduce(function(x,y) {merge(x,y)}, datalist)

########################################
data1 <- MeanExpForClusters1
data <- MeanExpForClusters
rownames(data) <- data$gene
data <- data[,-1]

rownames(data1) <- data1$gene
data1 <- data1[,-1]

######################################### write data both perecnt and mean
if (write.data == TRUE) {
MyDataToWrite <- cbind(data, data1)
MyDataToWrite <- cbind(rows = rownames(MyDataToWrite), MyDataToWrite)
  write.table((MyDataToWrite),file="bubble.gg.plot.data.tsv", sep="\t", row.names = FALSE)
}
if (write.data == FALSE) {
######################################### sacle
data <- as.matrix(data)
# data <- scale(log2(data + 1)) # this makes better cols but fails with moothing
data <- t(data)
data <- scale(data)
data <- t(data)
head(data)
# fix scale
FixScale <- function (mydata, min, max){
  Mydat <- mydata
  Mydat[Mydat > max] <- max
  Mydat[Mydat < min] <- min
  return(Mydat)
}
#
data <- FixScale(mydata = data, min = min.scale, max = max.scale)
########################################## merge Expression data and perecnt data
dat <- melt(data) # epression data
data1 <- as.matrix(data1) # perecnt data
dat1 <- melt(data1)
###
dat <- cbind(dat, Percent.Expressed = dat1$value)
colnames(dat) <- c("gene","cluster","Expression","Percent.Expressed")
# head(dat)
######################################### order genes and cluster order
dat$gene <- factor(dat$gene, gene)
cl <- as.character(unique(dat$cluster))
cl <- cl[order(nchar(cl))]
dat$cluster <- factor(dat$cluster, cl)
#####################################
# gene order
dat$gene <- factor(dat$gene,levels=rev(levels(dat$gene)))

############################################ plot
heatmap <- ggplot(dat, aes(x =  cluster, y = gene,
  colour = get(colour), size = get(size))) +
  geom_point() +
  scale_size(size) +
  scale_colour_gradient(low = heat.colors[1], high = heat.colors[2], name=colour) +
  theme_bw() + theme(axis.text.x=element_text(angle=90))
###
if (interactive == TRUE) {
  OUT.PUT <- paste(out.name, ".html", sep="")
  htmlwidgets::saveWidget(ggplotly(heatmap), OUT.PUT)
} else {
  return(heatmap)
  }
 }
}






