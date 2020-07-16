#' Impute data
#'
#' This function imputes data.
#' @param x An object of class iCellR.
#' @return An object of class iCellR
#' @export
myImp <- function (x = NULL) {
  if ("iCellR" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  ### get ADT data
    DATA <- x@main.data
    DATA <- DATA[ rowSums(DATA) > 0, ]
    DATA <- as.data.frame(t(DATA))
    data <- as.data.frame(t(x@pca.data[1:5]))
    dists = as.data.frame(as.matrix(dist(t(data), method = "euclidean")))[1]
    colnames(dists) <- "myDist"
    mydata <- as.data.frame(cbind(dists,DATA))
    mydata <- (mydata[order(mydata$myDist, decreasing = TRUE),])
    mydata <-  mydata[,-1]
    LG <- length(row.names(mydata))
    PRCENT = round((LG/100) * 4)
    TIMES <- round(LG / PRCENT)
    KK = sort(rep(1:PRCENT,TIMES))
    KK2 = rep(max(KK) + 1,LG - length(KK))
    KK3 = c(KK,KK2)
    df <- aggregate(x = mydata, by = list(KK3), FUN = "mean")
    df <- df[rep(seq_len(nrow(df)), each=TIMES),]
    df <- as.data.frame(df)
    df <- head(df,length(row.names(mydata)))
    df <-  df[,-1]
    row.names(df) <- row.names(mydata)
    ORDER <- as.data.frame(cbind(row.names(DATA),c(1:LG)))
    row.names(ORDER) <- ORDER$V1
#    ha = merge(ORDER,df, by = row.names)
    dd <- as.data.frame(merge(ORDER,df,by="row.names"))
    dd <- dd[,-1]
    row.names(dd) <- dd$V1
    dd <- dd[,-1]
    dd <- dd[,-1]
    dd <- as.data.frame(t(dd))
#    aggregate(x = testDF, by = list(by1), FUN = "sum")
  ####
  attributes(x)$imputed.data <- dd
  return(x)
}
