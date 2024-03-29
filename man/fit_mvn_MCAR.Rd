% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mvn_MCAR.R
\name{fit_mvn_MCAR}
\alias{fit_mvn_MCAR}
\title{Multivariate normal spatial mixture model clustering}
\usage{
fit_mvn_MCAR(Y, coords_df, K, nsim = 2000, burn = 1000, z_init = NULL)
}
\arguments{
\item{Y}{An n x g matrix of gene expression values. n is the number of cell spots and g is the number of features.}

\item{coords_df}{An n x 2 data frame or matrix of 2d spot coordinates.}

\item{K}{The number of mixture components to fit.}

\item{nsim}{Number of total MCMC iterations to run.}

\item{burn}{Number of MCMC iterations to discard as burn in. The number of saved samples is nsim - burn.}

\item{z_init}{Optional initialized allocation vector. Randomly initialized if NULL.}
}
\value{
a list of posterior samples
}
\description{
Implement Gibbs sampling for MVN model with MCAR spatial random effects
}
\examples{
\dontrun{
# parameters
data(coords_df_sim)
coords_df <- coords_df_sim[,1:2]
z <- remap_canonical2(coords_df_sim$z)
W_nn <- scran::buildKNNGraph(as.matrix(coords_df),
                             k = 4,
                             transposed = TRUE)
A <- igraph::as_adjacency_matrix(W_nn,
                                 type = "both",
                                 sparse = FALSE)
                                 
n <- nrow(coords_df) # number of observations
g <- 2 # number of features
K <- length(unique(coords_df_sim$z)) # number of clusters (mixture components)
pi <- table(z)/length(z) # cluster membership probability

# Cluster Specific Parameters
Mu <- list(
  Mu1 = rnorm(g,-5,1),
  Mu2 = rnorm(g,0,1),
  Mu3 = rnorm(g,5,1),
  Mu4 = rnorm(g,-2,3)
)
# cluster specific variance-covariance
S <- matrix(1,nrow = g,ncol = g) # y covariance matrix
diag(S) <- 1.5
Sig <- list(
  Sig1 = S,
  Sig2 = S, 
  Sig3 = S,
  Sig4 = S
)

# generate phi - not cluster specific
# conditional covariance of phi_i given phi_noti
m <- colSums(A)
M <- diag(m)
V <- matrix(0.4,nrow = g, ncol = g) # CAR covariance
diag(V) <- 0.6
V_true <- V
rho <- 0.999999 # Spatial dependence parameter ~ 1 for intrinsic CAR
Q <- diag(m) - rho*A	# m is number of neighbors for each spot
covphi <- solve(Q) \%x\% V # gn x gn covariance of phis
phi <- mvtnorm::rmvnorm(1, sigma=covphi)		# gn vector of spatial effects
PHI <- matrix(phi, ncol=g, byrow=TRUE) 	# n x g matrix of spatial effects
PHI <- t(scale(t(PHI)))

Y <- matrix(0, nrow = n, ncol = g)
for(i in 1:n)
{
  Y[i,] <- mvtnorm::rmvnorm(1,mean = Mu[[z[i]]] + PHI[i,],sigma = Sig[[z[i]]])
}

# fit model
# in practice use more mcmc iterations
fit_MCAR <- fit_mvn_MCAR(Y = Y, coords_df = coords_df, K = K, nsim = 10, burn = 0)}
}
