setClass("iCellR", representation (raw.data = "data.frame",
                                   main.data = "data.frame",
                                   imputed.data = "data.frame",
                                   scaled.data = "data.frame",
                                   batch.aligned.data = "data.frame",
                                   metadata = "data.frame",
                                   atac.raw = "data.frame",
                                   atac.main = "data.frame",
                                   spatial.data = "data.frame",
                                   stats = "data.frame",
                                   my.filters = "character",
                                   gene.model = "character",
                                   my.freq = "data.frame",
                                   obj.info = "character",
                                   vdj.data = "data.frame",
                                   tsne.data = "data.frame",
                                   tsne.data.3d = "data.frame",
                                   knetl.data = "data.frame",
                                   knetl.data.3d = "data.frame",
                                   umap.data = "data.frame",
                                   diffusion.data = "data.frame",
                                   pca.data = "data.frame",
                                   pca.info = "list",
                                   cca.data = "data.frame",
                                   opt.pcs = "numeric",
                                   dist.data = "data.frame",
                                   pseudo.mapA = "data.frame",
                                   pseudo.mapB = "data.frame",
                                   clust.avg = "data.frame",
                                   gene.data = "data.frame",
                                   adt.raw = "data.frame",
                                   adt.main = "data.frame",
                                   clust.cond.freq = "data.frame",
                                   cluster.data = "list",
                                   best.clust = "data.frame",
                                   data.conditions = "character",
                                   norm.factors = "data.frame",
                                   extra.data1 = "data.frame",
                                   extra.data2 = "data.frame"))
# hide slots
setMethod("show",
          "iCellR",
          function(object){
            message("###################################","")
            message(",--. ,-----.       ,--.,--.,------. ","")
            message("`--''  .--./ ,---. |  ||  ||  .--. ' ","")
            message(",--.|  |    | .-. :|  ||  ||  '--'.' ","")
            message("|  |'  '--'\\   --. |  ||  ||  |\  \ ","")
            message("`--' `-----' `----'`--'`--'`--' '--' ","")
            message("###################################","")
            message(object@obj.info,"")
            message("###################################","")
            message("   QC stats performed:",dim(object@stats)[1] != 0,", ","PCA performed:",dim(object@pca.data)[1] != 0)
            message("   Clustering performed:",!is.null(object@best.clust$clusters),", ","Number of clusters:",length(unique(object@best.clust$clusters)))
            message("   tSNE performed:",dim(object@tsne.data)[1] != 0,", ","UMAP performed:",dim(object@umap.data)[1] != 0,", ","DiffMap performed:",dim(object@diffusion.data)[1] != 0)
            message("   Main data dimensions (rows,columns): ",dim(object@main.data)[1],",",dim(object@main.data)[2])
            if (dim(object@main.data)[1] != 0){
              Cells <- colnames(object@main.data)
              InFO <- as.character((unique(data.frame(do.call('rbind', strsplit(as.character(Cells),'_',fixed=TRUE)))[1]))$X1)
              if(length(InFO) != 0) {
                InFO <- table(data.frame(do.call('rbind', strsplit(as.character(colnames(object@main.data)),'_',fixed=TRUE)))[1])
                D1 <- paste(as.character(as.data.frame(InFO)$Var1), collapse=",")
                D2 <- paste(as.character(as.data.frame(InFO)$Freq), collapse=",")
                message("   Data conditions in main data:",D1,"(", D2,")")
              }
            }
            message("   Normalization factors:",head(object@norm.factors,1),",","... ")
            message("   Imputed data dimensions (rows,columns):",dim(object@imputed.data)[1],",",dim(object@imputed.data)[2])
            message("############## scVDJ-seq ###########","")
            message("VDJ data dimentions (rows,columns):",dim(object@vdj.data)[1],",",dim(object@vdj.data)[2])
            message("############## CITE-seq ############","")
            message("   ADT raw data  dimensions (rows,columns):",dim(object@adt.raw)[1],",",dim(object@adt.raw)[2])
            message("   ADT main data  dimensions (rows,columns):",dim(object@adt.main)[1],",",dim(object@adt.main)[2])
            message("   ADT columns names:",head(colnames(object@adt.main),1),"... ")
            message("   ADT row names:",head(row.names(object@adt.main),1),"... ")
            message("############## scATAC-seq ############","")
            message("   ATAC raw data  dimensions (rows,columns):",dim(object@atac.raw)[1],",",dim(object@atac.raw)[2])
            message("   ATAC main data  dimensions (rows,columns):",dim(object@atac.main)[1],",",dim(object@atac.main)[2])
            message("   ATAC columns names:",head(colnames(object@atac.main),1),"... ")
            message("   ATAC row names:",head(row.names(object@atac.main),1),"... ")
            message("############## Spatial ###########","")
            message("Spatial data dimentions (rows,columns):",dim(object@spatial.data)[1],",",dim(object@spatial.data)[2])
            message("########### iCellR object ##########","")
          })
