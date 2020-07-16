#' K Nearest Neighbour Search
#'
#' Uses a kd-tree to find the p number of near neighbours for each point in an input/output dataset.
#'
#' Use the nn2 function from the RANN package, utilizes the Approximate Near Neighbor (ANN) C++ library,
#' which can give the exact near neighbours or (as the name suggests) approximate near neighbours
#' to within a specified error bound. For more information on the ANN library please
#' visit http://www.cs.umd.edu/~mount/ANN/.
#'
#' @param data matrix; input data matrix
#' @param k integer; number of nearest neighbours
#' @importFrom RANN nn2
#' @export
find_neighbors <- function(data, k){
  nearest <- nn2(data, data, k, searchtype = "standard")
  return(nearest[[1]])
}
