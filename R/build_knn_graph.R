#' Make KNN network
#'
#' Construct a binary adjacency matrix
#'
#' @param coords_df An n x 2 data frame or matrix of 2d spot coordinates  
#' @param k The number of neighbors
#'
#' @return an adjacency matrix
#' @export
#' @importFrom igraph graph_from_edgelist as_adjacency_matrix
#' @examples 
#' data(coords_df_sim)
#' coords_df <- coords_df_sim[,1:2]
build_knn_graph <- function(coords,k)
{
  dist <- as.matrix(dist(as.matrix(coords)))
  edges <- mat.or.vec(0,2)
  
  for (i in 1:nrow(dist))
  {
    matches <- setdiff(order(dist[i,],decreasing = FALSE)[1:(k+1)],i)
    edges <- rbind(edges,cbind(rep(i,k),matches))  
    edges <- rbind(edges,cbind(matches,rep(i,k)))  
  }
  
  G <- igraph::graph_from_edgelist(edges,directed=FALSE)
  A <- igraph::as_adjacency_matrix(G,type = "both",sparse = FALSE)
  A[A > 0] <- 1
  return(A)        
}