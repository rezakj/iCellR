#' Calculate Cell cycle phase prediction
#'
#' This function takes an object of class iCellR and assignes cell cycle stage for the cells.
#' @param object A data frame containing gene counts for cells.
#' @param s.genes Genes that are used as a marker for S phase.
#' @param g2m.genes Genes that are used as a marker for G2 and M phase.
#' @return The data frame object
#' @importFrom Hmisc cut2
#' @export
cc <- function (object = NULL,
                s.genes = s.phase,
                g2m.genes = g2m.phase) {
  if ("iCellR" != class(object)[1]) {
    stop("object should be an object of class iCellR")
  }
  ##### get genes case insensetive
  ALLgenes = row.names(object@raw.data)
  s.phase.genes <- s.genes
  s.phase.genes <- paste("^",s.phase.genes,"$", sep="")
  s.phase.genes <- paste(s.phase.genes,collapse="|")
  s.phase.genes <- grep(s.phase.genes, x = ALLgenes, value = TRUE, ignore.case = TRUE)
  s.genes <- s.phase.genes
  #
  g2m.phase.genes <- g2m.genes
  g2m.phase.genes <- paste("^",g2m.phase.genes,"$", sep="")
  g2m.phase.genes <- paste(g2m.phase.genes,collapse="|")
  g2m.phase.genes <- grep(g2m.phase.genes, x = ALLgenes, value = TRUE, ignore.case = TRUE)
  g2m.genes <- g2m.phase.genes
  #####
  Table <-  object@stats
  row.names(Table) <- Table$CellIds
  attributes(object)$stats <- Table
#  head(object@stats)
  ########## 1
  LengthCheck <- function(values, cutoff = 0) {
    return(vapply(
      X = values,
      FUN = function(x) {
        return(length(x = x) > cutoff)
      },
      FUN.VALUE = logical(1)
    ))
  }

  ########## 2
    ####
    AddMetaDatame <- function (object, metadata, col.name = NULL)
    {
      if (typeof(x = metadata) != "list") {
        metadata <- as.data.frame(x = metadata)
        if (is.null(x = col.name)) {
          stop("Please provide a name for provided metadata")
        }
        colnames(x = metadata) <- col.name
      }
      cols.add <- colnames(x = metadata)
      meta.order <- match(rownames(object@stats), rownames(metadata))
      meta.add <- metadata[meta.order, ]
      if (all(is.null(x = meta.add))) {
        stop("Metadata provided doesn't match the cells in this object")
      }
      object@stats[, cols.add] <- meta.add
      return(object)
    }
    ######### func 2
    AddModuleScoreme <- function (object, genes.list = NULL, genes.pool = NULL, n.bin = 25,
              seed.use = 1, ctrl.size = 100, use.k = FALSE, enrich.name = "Cluster",
              random.seed = 1)
    {
      set.seed(seed = random.seed)
      genes.old <- genes.list
      if (use.k) {
        genes.list <- list()
        for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
          genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster ==
                                               i))
        }
        cluster.length <- length(x = genes.list)
      }
      else {
        if (is.null(x = genes.list)) {
          stop("Missing input gene list")
        }
        genes.list <- lapply(X = genes.list, FUN = function(x) {
          return(intersect(x = x, y = rownames(x = object@raw.data)))
        })
        cluster.length <- length(x = genes.list)
      }

      if (!all(LengthCheck(values = genes.list))) {
        warning(paste("Could not find enough genes in the object from the following gene lists:",
                      paste(names(x = which(x = !LengthCheck(values = genes.list)))),
                      "Attempting to match case..."))
        genes.list <- lapply(X = genes.old, FUN = CaseMatch,
                             match = rownames(x = object@raw.data))
      }
      if (!all(LengthCheck(values = genes.list))) {
        stop(paste("The following gene lists do not have enough genes present in the object:",
                   paste(names(x = which(x = !LengthCheck(values = genes.list)))),
                   "exiting..."))
      }

      if (is.null(x = genes.pool)) {
        genes.pool = rownames(x = object@raw.data)
      }
      data.avg <- Matrix::rowMeans(x = object@raw.data[genes.pool,
                                                   ])
      data.avg <- data.avg[order(data.avg)]
      #
      data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg)/n.bin)))
      names(x = data.cut) <- names(x = data.avg)
      ctrl.use <- vector(mode = "list", length = cluster.length)
      for (i in 1:cluster.length) {
        genes.use <- genes.list[[i]]
        for (j in 1:length(x = genes.use)) {
          ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut ==
                                                                                  data.cut[genes.use[j]])], size = ctrl.size, replace = FALSE)))
        }
      }
      ctrl.use <- lapply(X = ctrl.use, FUN = unique)
      ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use),
                            ncol = ncol(x = object@raw.data))
      for (i in 1:length(ctrl.use)) {
        genes.use <- ctrl.use[[i]]
        ctrl.scores[i, ] <- Matrix::colMeans(x = object@raw.data[genes.use,
                                                             ])
      }
      genes.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length,
                             ncol = ncol(x = object@raw.data))
      for (i in 1:cluster.length) {
        genes.use <- genes.list[[i]]
        data.use <- object@raw.data[genes.use, , drop = FALSE]
        genes.scores[i, ] <- Matrix::colMeans(x = data.use)
      }
      genes.scores.use <- genes.scores - ctrl.scores
      rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
      genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))
      rownames(x = genes.scores.use) <- colnames(x = object@raw.data)
      object <- AddMetaDatame(object = object, metadata = genes.scores.use,
                            col.name = colnames(x = genes.scores.use))
      gc(verbose = FALSE)
      return(object)
    }
    ######## 3
    CellCycleScoringme <- function (object, g2m.genes, s.genes, set.ident = FALSE)
    {
      enrich.name <- "Cell Cycle"
      genes.list <- list(S.Score = s.genes, G2M.Score = g2m.genes)
      object.cc <- AddModuleScoreme(object = object, genes.list = genes.list,
                                    enrich.name = enrich.name, ctrl.size = min(vapply(X = genes.list,
                                                                                      FUN = length, FUN.VALUE = numeric(1))))
      cc.columns <- grep(pattern = enrich.name, x = colnames(x = object.cc@stats))
      cc.scores <- object.cc@stats[, cc.columns]
#      head(cc.scores)
      rm(object.cc)
      gc(verbose = FALSE)
      assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores,
                                                                     first = "S", second = "G2M", null = "G1") {
        if (all(scores < 0)) {
          return(null)
        }
        else {
          if (length(which(x = scores == max(scores))) > 1) {
            return("Undecided")
          }
          else {
            return(c(first, second)[which(x = scores == max(scores))])
          }
        }
      })
      cc.scores <- merge(x = cc.scores, y = data.frame(assignments),
                         by = 0)
      colnames(x = cc.scores) <- c("rownames", "S.Score", "G2M.Score",
                                   "Phase")
      rownames(x = cc.scores) <- cc.scores$rownames
      cc.scores <- cc.scores[, c("S.Score", "G2M.Score", "Phase")]
      #    cc.scores
      object <- AddMetaDatame(object = object, metadata = cc.scores)
      if (set.ident) {
        object <- StashIdent(object = object, save.name = "old.ident")
        object <- SetAllIdent(object = object, id = "Phase")
      }
      return(object)
}
#############
#############
############# How to run
    object <- CellCycleScoringme(object = object, s.genes = s.genes, g2m.genes = g2m.genes)
    STATS <- object@stats
  attributes(object)$stats <- STATS
  return(object)
}
