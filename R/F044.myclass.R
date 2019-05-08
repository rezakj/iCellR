setClass("iCellR", representation (raw.data = "data.frame",
                                   main.data = "data.frame",
                                   imputed.data = "data.frame",
                                   scaled.data = "data.frame",
                                   stats = "data.frame",
                                   obj.info = "character",
                                   vdj.data = "data.frame",
                                   tsne.data = "data.frame",
                                   umap.data = "data.frame",
                                   diffusion.data = "data.frame",
                                   pca.data = "data.frame",
                                   pca.info = "list",
                                   cca.data = "data.frame",
                                   opt.pcs = "numeric",
                                   dist.data = "data.frame",
                                   diff.st.data = "data.frame",
                                   tsne.data.3d = "data.frame",
                                   clust.avg = "data.frame",
                                   gene.data = "data.frame",
                                   gene.model = "data.frame",
                                   adt.raw = "data.frame",
                                   adt.main = "data.frame",
                                   clust.cond.freq = "data.frame",
                                   cluster.data = "list",
                                   best.clust = "data.frame",
                                   data.conditions = "character",
                                   norm.factors = "data.frame"))
# hide slots
setMethod("show",
          "iCellR",
          function(object){
            cat("###################################","\n")
            cat(",--. ,-----.       ,--.,--.,------. ","\n")
            cat("`--''  .--./ ,---. |  ||  ||  .--. ' ","\n")
            cat(",--.|  |    | .-. :|  ||  ||  '--'.' ","\n")
            cat("|  |'  '--'\\   --. |  ||  ||  |\  \ ","\n")
            cat("`--' `-----' `----'`--'`--'`--' '--' ","\n")
            cat("###################################","\n")
            cat(object@obj.info[1],"\n")
            cat("   ",object@obj.info[2],"\n")
            cat("   ",object@obj.info[3],"\n")
            cat("   ",object@obj.info[4],"\n")
            cat("   ",object@obj.info[5],"\n")
            cat("###################################","\n")
            cat("   QC stats performed:",!is.null(object@stats),", ")
            cat("PCA performed:",!is.null(object@pca.info),", ")
            cat("CCA performed:",!is.null(object@cca.data),"\n")
            cat("   Clustering performed:",!is.null(object@best.clust$clusters),", ")
            cat("Number of clusters:",length(unique(object@best.clust$clusters)),"\n")
            cat("   tSNE performed:",!is.null(object@tsne.data),", ")
            cat("UMAP performed:",!is.null(object@umap.data),", ")
            cat("DiffMap performed:",!is.null(object@diffusion.data),"\n")
            cat("   Main data dimentions (rows,columns):",dim(object@main.data),"\n")
            cat("   Normalization factors:",head(object@norm.factors,3),"... \n")
            cat("   Imputed data dimentions (rows,columns):",dim(object@imputed.data),"\n")
            cat("############## scVDJ-Seq ###########","\n")
            cat("   VDJ data dimentions (rows,columns):",dim(object@vdj.data),"\n")
            cat("############## CITE-Seq ############","\n")
            cat("   ADT raw data dimentions (rows,columns):",dim(object@adt.raw),"\n")
            cat("   ADT main data dimentions (rows,columns):",dim(object@adt.main),"\n")
            cat("   ADT columns names:",head(colnames(object@adt.main),3),"... \n")
            cat("   ADT row names:",head(row.names(object@adt.main),3),"... \n")
            cat("######## iCellR object made ########","\n")
          })
