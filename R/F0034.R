#' Differential expression (DE) analysis
#'
#' This function takes an object of class iCellR and performs differential expression (DE) analysis for clusters and conditions.
#' @param x An object of class iCellR.
#' @param data.type Choose from "main" and "imputed", default = "main"
#' @param de.by Choose from "clusters", "conditions", "clustBase.condComp" or "condBase.clustComp".
#' @param cond.1 First condition to do DE analysis on.
#' @param cond.2 Second condition to do DE analysis on.
#' @param base.cond A base condition or cluster if de.by is either cond.clust or clust.cond
#' @return An object of class iCellR
#' @examples
#' diff.res <- run.diff.exp(demo.obj, de.by = "clusters", cond.1 = c(1), cond.2 = c(2))
#'
#' head(diff.res)
#'
#' @export
run.diff.exp <- function (x = NULL,
                      data.type = "main",
                      de.by = "clusters",
                      cond.1 = "array",
                      cond.2 = "array",
                      base.cond = 0) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ###########
#  dat <- x@main.data
  ## get main data
  if (data.type == "main") {
    dat <- x@main.data
  }
  if (data.type == "imputed") {
    dat <- x@imputed.data
  }
  # 2 dimentions
  DATA <- x@best.clust
  ############## set wich clusters you want as condition 1 and 2
  CondA = cond.1
  CondB = cond.2
  CondAnames = paste(CondA, collapse="_")
  CondBnames = paste(CondB, collapse="_")
  if (de.by == "clustBase.condComp") {
    CondAnames = paste(CondA, collapse="_")
    CondAnames = paste(CondAnames,"inClust",base.cond,sep=".")
    CondBnames = paste(CondB, collapse="_")
    CondBnames = paste(CondBnames,"inClust",base.cond,sep=".")
  }
  if (de.by == "condBase.clustComp") {
    CondAnames = paste(CondA, collapse="_")
    CondAnames = paste(CondAnames,"inCond",base.cond,sep=".")
    CondBnames = paste(CondB, collapse="_")
    CondBnames = paste(CondBnames,"inCond",base.cond,sep=".")
  }
  # conditions
  if (de.by == "conditions") {
    col.legend <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
    conditions <- as.character(as.matrix(col.legend))
    Table <- as.data.frame(cbind(DATA, conditions = conditions))
    Cluster0 <- row.names(subset(Table, Table$conditions %in% CondA))
    Cluster1 <- row.names(subset(Table, Table$conditions %in% CondB))
  }
  # clusters
  if (de.by == "clusters") {
    if (is.null(DATA$clusters)) {
      stop("Clusters are not assigend yet, please run assign.clust fisrt.")
    } else {
      Table=DATA
      Cluster0 <- row.names(subset(Table, Table$clusters %in% CondA))
      Cluster1 <- row.names(subset(Table, Table$clusters %in% CondB))
    }
  }
  #
  if (de.by == "clustBase.condComp") {
    if (base.cond == 0) {
      stop("You should choose a base cluster/clusters to compare the conditions in it/them")
    }
    my.base.cond = base.cond
    col.legend <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
    conditions <- as.character(as.matrix(col.legend))
    Table <- as.data.frame(cbind(DATA, conditions = conditions))
    Table <- subset(Table, Table$clusters %in% my.base.cond)
    Cluster0 <- row.names(subset(Table, Table$conditions %in% CondA))
    Cluster1 <- row.names(subset(Table, Table$conditions %in% CondB))
  }
  #
  if (de.by == "condBase.clustComp") {
    if (base.cond == 0) {
      stop("You should choose a base condition/conditions to compare the clusters in it/them")
    }
    my.base.cond = base.cond
    col.legend <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
    conditions <- as.character(as.matrix(col.legend))
    Table <- as.data.frame(cbind(DATA, conditions = conditions))
    Table <- subset(Table, Table$conditions %in% my.base.cond)
    Cluster0 <- row.names(subset(Table, Table$clusters %in% CondA))
    Cluster1 <- row.names(subset(Table, Table$clusters %in% CondB))
  }
  ############## Filter
  cond1 <- dat[,Cluster0]
  cond2 <- dat[,Cluster1]
  ### merge both for pval length not matching error
  mrgd <- cbind(cond1,cond2)
  #    mrgd <- merge(cond1, cond2, by="row.names")
  #    row.names(mrgd) <- mrgd$Row.names
  #    mrgd <- mrgd[,-1]
  mrgd <- data.matrix(mrgd)
  # mean
  meansCond1 <- apply(cond1, 1, function(cond1) {mean(cond1)})
  meansCond2 <- apply(cond2, 1, function(cond2) {mean(cond2)})
  baseMean <- apply(mrgd, 1, function(mrgd) {mean(mrgd)})
  # FC
  FC <- meansCond2/meansCond1
  FC.log2 <- log2(FC)
  # dims
  Cond1_Start <- 1
  Cond1_End <- dim(cond1)[2]
  Cond2_Start <- dim(cond1)[2] + 1
  Cond2_End <- dim(cond1)[2] + dim(cond2)[2]
  # pval
  Pval <- apply(mrgd, 1, function(mrgd) {
    t.test(x = mrgd[Cond1_Start:Cond1_End], y = mrgd[Cond2_Start:Cond2_End])$p.value
  })
  # padj
  FDR <- p.adjust(Pval)
  # combine
  Stats <- cbind(
    baseMean = baseMean,
    MeanCond1 = meansCond1,
    MeanCond2 = meansCond2,
    FC=FC,
    FC.log2=FC.log2,
    t.test=Pval,
    t.test.adj=FDR)
  # name column
  colnames(Stats) <- c("baseMean",CondAnames, CondBnames,"foldChange","log2FoldChange","pval","padj")
  # return
  return(Stats)
}
