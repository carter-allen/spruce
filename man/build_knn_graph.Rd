% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_knn_graph.R
\name{build_knn_graph}
\alias{build_knn_graph}
\title{Make KNN network}
\usage{
build_knn_graph(coords, k)
}
\arguments{
\item{coords}{An n x 2 data frame or matrix of 2d spot coordinates}

\item{k}{The number of neighbors}
}
\value{
an adjacency matrix
}
\description{
Construct a binary adjacency matrix
}
\examples{
data(coords_df_sim)
coords_df <- coords_df_sim[,1:2]
A <- build_knn_graph(coords_df,k = 4)
}
