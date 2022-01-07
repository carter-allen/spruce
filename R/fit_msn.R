#' Multivariate skew-normal mixture model clustering
#'
#' Implement Gibbs sampling for MSN model with no spatial random effects
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
#' @importFrom stats cov
#' @importFrom MCMCpack rdirichlet
#' @importFrom truncnorm rtruncnorm
#' @importFrom Rclusterpp Rclusterpp.hclust
#' @examples
#' \dontrun{ 
#' # parameters
#' n <- 1000 # 
#' number of observations
#' g <- 3 # number of features
#' K <- 3 # number of clusters (mixture components)
#' pi <- rep(1/K,K) # cluster membership probability
#' z <- sample(1:K, size = n, replace = TRUE, prob = pi) # cluster indicators
#' z <- remap_canonical2(z)
#' t_true <- truncnorm::rtruncnorm(n,0,Inf,0,1)
#' t <- t_true
#' 
#' # Cluster Specific Parameters
#' # cluster specific means
#' Mu <- list(
#'   Mu1 = rnorm(g,-5,1),
#'   Mu2 = rnorm(g,0,1),
#'   Mu3 = rnorm(g,5,1)
#' )
#' # Cluster speficic skewness
#' Xi <- list(
#'   Xi1 = rep(2,g),
#'   Xi2 = rep(0,g),
#'   Xi3 = rep(-3,g)
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
#'   Y[i,] <- mvtnorm::rmvnorm(1,mean = Mu[[z[i]]] + t[i]*Xi[[z[i]]],sigma = Sig[[z[i]]])
#' }
#' 
#' # fit model
#' fit1 <- fit_msn_clustering(Y,3,10,0)}

fit_msn <- function(Y,
                    K,
                    nsim = 2000,
                    burn = 1000,
                    z_init = NULL)
{
  # parameters
  n <- nrow(Y) # number of observations
  p <- ncol(Y) # number of features
  pi <- rep(1/K,K) # cluster membership probability
  if(is.null(z_init)) # initialize z
  {
    fit_hclust <- Rclusterpp::Rclusterpp.hclust(Y)
    z_init <- stats::cutree(fit_hclust,k = K)
    z <- z_init
  }
  else # user provided initialization
  {
    z <- z_init
    pi <- table(z)/n
  }
  ts <- truncnorm::rtruncnorm(n,0,Inf,0,1)
  
  # priors - shared across clusters
  mu0 <- rep(0,p)
  xi0 <- rep(0,p)
  L0 <- S0 <- P0 <- diag(p)
  nu0 <- 2
  a0 <- rep(4,K) # prior parameter vector for pi1,...,piK
  
  # cluster specific sample stats
  Sigma <- list(0)
  Ybar <- list(0)
  for(k in 1:K)
  {
    Sigma[[k]] <- stats::cov(Y[z == k,])
    Ybar[[k]] <- colMeans(Y[z == k,])
  }
  
  # Intermediate MCMC vars
  Ln <- Pn <- list(0)
  mn <- xn <- list(0)
  mun <- list(0)
  xin <- list(0)
  Sn <- list(0)
  
  # Empty sample storage
  MU <- XI <- SIGMA <- vector("list",K)
  
  n_save <- nsim - burn
  Z <- matrix(0,nrow = n_save,ncol = n)
  for(k in 1:K)
  {
    mun[[k]] <- rep(0,p)
    xin[[k]] <- rep(0,p)
    MU[[k]] <- matrix(0,nrow = n_save,ncol = p)
    XI[[k]] <- matrix(0,nrow = n_save,ncol = p)
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
      tk <- ts[z == k]
      Yk <- Y[z == k,]
      tkmat <- matrix(tk,nrow = nk, ncol = 1)
      Etk <- Yk - tkmat %*% xin[[k]]
      Etk_bar <- colMeans(Etk)
      
      ### update mu - cluster specific
      Ln[[k]] <- solve(solve(L0) + nk*solve(Sigma[[k]]))
      mn[[k]] <- Ln[[k]] %*% (solve(L0) %*% mu0 + nk*solve(Sigma[[k]]) %*% Etk_bar) 
      mun[[k]] <- mvrnormArma(1,mn[[k]],Ln[[k]])
      Munk <- matrix(mun[[k]],nrow = nk,ncol = p,byrow = TRUE)
      
      ### update xi - cluster specific
      Pn[[k]] <- solve(solve(P0) + sum(tk^2)*solve(Sigma[[k]]))
      xn[[k]] <- Pn[[k]] %*% (solve(P0) %*% xi0 + solve(Sigma[[k]]) %*% colSums(tk * (Yk - Munk)))
      xin[[k]] <- mvrnormArma(1,xn[[k]],Pn[[k]])
      
      ### update Sigma - cluster specific 
      Ek <- Etk - Munk
      Sn[[k]] <- S0 + t(Ek) %*% Ek 
      Sigma[[k]] <- solve(r2arma::rwishart(nu0+nk, solve(Sn[[k]])))
      
      ### update t
      k_inds <- (1:n)[z == k]
      Ak <- solve(1 + t(xn[[k]]) %*% solve(Sigma[[k]]) %*% xn[[k]])
      ts <- update_t(ts,k,k_inds,Ak,xn[[k]],Sigma[[k]],Y,mun[[k]],0,Inf)
    }
    
    ### Update cluster indicators
    z <- update_z_MSN(z,Y,ts,mun,xin,Sigma,pi,1:K)
    # remap to address label switching
    z <- remap_canonical2(z)
    n.z <- as.vector(unname(table(z))) # gives the number of members currently in each class
    
    # Update pi1,...,piK 
    pi <- MCMCpack::rdirichlet(1,a0 + n.z)
    
    ## save results
    if(i > burn)
    {
      iter <- i - burn
      for(k in 1:K)
      {
        MU[[k]][iter,] <- mun[[k]]
        XI[[k]][iter,] <- xin[[k]]
        SIGMA[[k]][iter,] <- c(Sigma[[k]])
      }
      Z[iter,] <- z
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  run.time<-proc.time()-start.time
  print(paste("Finished MCMC after",run.time[1],"seconds"))
  
  z_map <- apply(Z, 2, get_map)
  
  ret_list <- list(Y = Y,
                   MU = MU,
                   XI = XI,
                   SIGMA = SIGMA,
                   K = K,
                   Z = Z,
                   z = z_map)
  return(ret_list)
}