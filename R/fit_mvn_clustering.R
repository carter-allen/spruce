#' Multivariate normal mixture model clustering
#'
#' Implement Gibbs sampling for MVN model with no spatial random effects
#'
#' @param Y An n x g matrix of gene expression values. n is the number of cell spots and g is the number of features.
#' @param K The number of mixture components to fit. 
#' @param nsim Number of total MCMC iterations to run.
#' @param burn Number of MCMC iterations to discard as burn in. The number of saved samples is nsim - burn.
#' @param z_init Optional initialized allocation vector. Randomly initialized if NULL. 
#'
#' @return a list of posterior samples
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom mvtnorm rmvnorm
#' @examples 
#' # parameters
#' n <- 1000 # number of observations
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
#' fit1 <- fit_mvn_clustering(Y,3,1000,100)

fit_mvn_clustering <- function(Y,K,nsim = 2000,burn = 1000,z_init = NULL)
{
  # parameters
  n <- nrow(Y) # number of observations
  p <- ncol(Y) # number of features
  pi <- rep(1/K,K) # cluster membership probability
  if(is.null(z_init))
  {
    z <- sample(1:K, size = n, replace = TRUE, prob = pi) # cluster indicators
    z <- remap_canonical2(z)
  }
  else
  {
    z <- z_init
  }
  
  # priors - shared across clusters
  mu0 <- colMeans(Y)
  L0 <- S0 <- diag(p)
  nu0 <- 2
  
  # cluster specific sample stats
  Sigma <- list(0)
  Ybar <- list(0)
  for(k in 1:K)
  {
    Sigma[[k]] <- cov(Y[z == k,])
    Ybar[[k]] <- colMeans(Y[z == k,])
  }
  
  # Intermediate MCMC vars
  Ln <- list(0)
  mn <- list(0)
  mun <- list(0)
  Sn <- list(0)
  
  # Empty sample storage
  MU <- SIGMA <- vector("list",K)
  
  n_save <- nsim - burn
  Z <- matrix(0,nrow = n_save,ncol = n)
  P <- matrix(0,nrow = n_save,ncol = n)
  for(k in 1:K)
  {
    MU[[k]] <- matrix(0,nrow = n_save,ncol = p)
    SIGMA[[k]] <- matrix(0,nrow = n_save,ncol = p*p)
  }
  start.time <- proc.time()
  print(paste("Started MCMC of",nsim))
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for(i in 1:nsim) 
  {
    ### Update cluster - specific parameters
    for(k in 1:K)
    {
      ### update cluster specific sample stats
      nk <- sum(z == k)
      Ybar[[k]] <- colMeans(Y[z == k,])
      
      ### update mu - cluster specific
      Ln[[k]] <- solve(solve(L0) + nk*solve(Sigma[[k]]))
      mn[[k]] <- Ln[[k]] %*% (solve(L0) %*% mu0 + nk*solve(Sigma[[k]]) %*% Ybar[[k]]) 
      mun[[k]] <- mvrnormArma(1 ,mn[[k]],Ln[[k]])
      
      ### update Sigma - cluster specific 
      Sn[[k]] <- S0 + (t(Y[z == k,]) - c(mun[[k]])) %*% t(t(Y[z == k,]) - c(mun[[k]])) 
      Sigma[[k]] <- solve(r2arma::rwishart(nu0+nk, solve(Sn[[k]])))
    }
    
    ### Update cluster indicators
    PZ <- update_z(z,Y,mun,Sigma,pi,1:K)
    z <- PZ[1,]
    # remap to address label switching
    z <- remap_canonical2(z)
    
    ## save results
    if(i > burn)
    {
      iter <- i - burn
      for(k in 1:K)
      {
        MU[[k]][iter,] <- mun[[k]]
        SIGMA[[k]][iter,] <- c(Sigma[[k]])
      }
      Z[iter,] <- z
      P[iter,] <- PZ[2,]
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  run.time<-proc.time()-start.time
  print(paste("Finished MCMC after",run.time[1],"seconds"))
  
  z_map <- apply(Z, 2, get_map)
  
  ret_list <- list(Y = Y,
                   MU = MU,
                   SIGMA = SIGMA,
                   K = K,
                   Z = Z,
                   z = z_map,
                   P = P)
  return(ret_list)
}