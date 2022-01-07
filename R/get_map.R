#' Get MAP estimate of cluster indicators
#'
#' Compute maximum aposteriori (MAP) estimate of cluster indicators
#'
#' @param z All cluster indicator posterior samples from a given cell spot
#'
#' @return MAP estimate of cluster labels. Useful applied over columns of posterior samples matrix (see example)
#' @export
#' @importFrom mvtnorm rmvnorm
#' @examples 
#' # parameters
#' n <- 100 # number of observations
#' g <- 3 # number of features
#' K <- 3 # number of clusters (mixture components)
#' pi <- rep(1/K,K) # cluster membership probability
#' z <- sample(1:K, size = n, replace = TRUE, prob = pi) # cluster indicators
#' z <- remap_canonical2(z)
#' 
#' # Cluster Specific Parameters
#' # cluster specific means
#' Mu <- list(
#'   Mu1 = rnorm(g,-5,1),
#'   Mu2 = rnorm(g,0,1),
#'   Mu3 = rnorm(g,5,1)
#' )
#' # cluster specific variance-covariance
#' S <- matrix(1,nrow = g,ncol = g) # covariance matrix
#' diag(S) <- 1.5
#' Sig <- list(
#'   Sig1 = S,
#'   Sig2 = S, 
#'   Sig3 = S
#' )
#' 
#' Y <- matrix(0, nrow = n, ncol = g)
#' for(i in 1:n)
#' {
#'   Y[i,] <- mvtnorm::rmvnorm(1,mean = Mu[[z[i]]],sigma = Sig[[z[i]]])
#' }
#' 
#' # fit model
#' fit1 <- fit_mvn(Y,3,100,0)
#'
#' # Apply get_map() to columns of Z (i.e., posterior samples from each cell spot)
#' z_map <- apply(fit1$Z, 2, get_map)

get_map <- function(z)
{
  return(names(which.max(table(z))))
}