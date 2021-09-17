#' Canonical re-mapping of mixture component labels
#'
#' Avoid label switching by re-mapping sampled mixture component labels at each iteration (Peng and Carvhalo 2016).
#'
#' @param z A length-n vector of discrete mixture component labels
#'
#' @return A length-n vector of mixture component labels re-mapped to a canonical sub-space
#' @export
#' @examples 
#' # parameters
#' n <- 10 # number of observations
#' K <- 3 # number of clusters (mixture components)
#' pi <- rep(1/K,K) # cluster membership probability
#' z <- sample(1:K, size = n, replace = TRUE, prob = pi) # cluster indicators
#' z <- remap_canonical2(z)

remap_canonical2 <- function(z)
{
  ord_obs <- unique(z)
  z_ret <- z
  for(k in 1:length(ord_obs))
  {
    z_ret[z == ord_obs[k]] <- k
  }
  return(z_ret)
}
