#' Multivariate normal spatial mixture model clustering
#'
#' Implement Gibbs sampling for MVN model with MCAR spatial random effects
#'
#' @param Y An n x g matrix of gene expression values. n is the number of cell spots and g is the number of features.
#' @param coords_df An n x 2 data frame or matrix of 2d spot coordinates.  
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
#' @importFrom scran buildKNNGraph
#' @importFrom igraph as_adjacency_matrix
#' @importFrom MCMCpack rdirichlet
#' @examples 
#' \dontrun{
#' # parameters
#' data(coords_df_sim)
#' coords_df <- coords_df_sim[,1:2]
#' z <- remap_canonical2(coords_df_sim$z)
#' W_nn <- scran::buildKNNGraph(as.matrix(coords_df),
#'                              k = 4,
#'                              transposed = TRUE)
#' A <- igraph::as_adjacency_matrix(W_nn,
#'                                  type = "both",
#'                                  sparse = FALSE)
#'                                  
#' n <- nrow(coords_df) # number of observations
#' g <- 2 # number of features
#' K <- length(unique(coords_df_sim$z)) # number of clusters (mixture components)
#' pi <- table(z)/length(z) # cluster membership probability
#' 
#' # Cluster Specific Parameters
# cluster specific means
#' Mu <- list(
#'   Mu1 = rnorm(g,-5,1),
#'   Mu2 = rnorm(g,0,1),
#'   Mu3 = rnorm(g,5,1),
#'   Mu4 = rnorm(g,-2,3)
#' )
#' # cluster specific variance-covariance
#' S <- matrix(1,nrow = g,ncol = g) # y covariance matrix
#' diag(S) <- 1.5
#' Sig <- list(
#'   Sig1 = S,
#'   Sig2 = S, 
#'   Sig3 = S,
#'   Sig4 = S
#' )
#' 
#' # generate phi - not cluster specific
#' # conditional covariance of phi_i given phi_noti
#' m <- colSums(A)
#' M <- diag(m)
#' V <- matrix(0.4,nrow = g, ncol = g) # CAR covariance
#' diag(V) <- 0.6
#' V_true <- V
#' rho <- 0.999999 # Spatial dependence parameter ~ 1 for intrinsic CAR
#' Q <- diag(m) - rho*A	# m is number of neighbors for each spot
#' covphi <- solve(Q) %x% V # gn x gn covariance of phis
#' phi <- mvtnorm::rmvnorm(1, sigma=covphi)		# gn vector of spatial effects
#' PHI <- matrix(phi, ncol=g, byrow=TRUE) 	# n x g matrix of spatial effects
#' PHI <- t(scale(t(PHI)))
#' 
#' Y <- matrix(0, nrow = n, ncol = g)
#' for(i in 1:n)
#' {
#'   Y[i,] <- mvtnorm::rmvnorm(1,mean = Mu[[z[i]]] + PHI[i,],sigma = Sig[[z[i]]])
#' }
#' 
#' # fit model
#' # in practice use more mcmc iterations
#' fit_MCAR <- fit_mvn_MCAR(Y = Y, coords_df = coords_df, K = K, nsim = 10, burn = 0)}

fit_mvn_MCAR <- function(Y,coords_df,K,nsim = 2000,burn = 1000,z_init = NULL)
{
  # parameters
  n <- nrow(coords_df) # number of observations
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
  
  # adjacency matrix
  W_nn <- scran::buildKNNGraph(as.matrix(coords_df),k = 4,transposed = TRUE)
  A <- igraph::as_adjacency_matrix(W_nn,type = "both",sparse = FALSE)
  m <- colSums(A)
  M <- diag(m)
  
  # random effects
  PHI <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  
  # priors - shared across clusters
  mu0 <- colMeans(Y)
  L0 <- S0 <- diag(p)
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
  Ln <- list(0)
  mn <- list(0)
  mun <- list(0)
  Sn <- list(0)
  
  # Empty sample storage
  MU <- SIGMA <- vector("list",K)
  
  n_save <- nsim - burn
  Z <- matrix(0,nrow = n_save,ncol = n)
  Phi <- PHI
  
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
      Ybar[[k]] <- colMeans(Y[z == k,] - Phi[z == k,])
      
      ### update mu - cluster specific
      Ln[[k]] <- solve(solve(L0) + nk*solve(Sigma[[k]]))
      mn[[k]] <- Ln[[k]] %*% (solve(L0) %*% mu0 + nk*solve(Sigma[[k]]) %*% Ybar[[k]]) 
      mun[[k]] <- mvrnormArma(1 ,mn[[k]],Ln[[k]])
      
      ### update Sigma - cluster specific 
      Sn[[k]] <- S0 + (t(Y[z == k,]) - t(Phi[z == k,]) - c(mun[[k]])) %*% t(t(Y[z == k,]) - t(Phi[z == k,]) - c(mun[[k]])) 
      Sigma[[k]] <- solve(r2arma::rwishart(nu0+nk, solve(Sn[[k]])))
    }
    
    ### Update random effects Phi
    Phi <- update_phi_spot_MCAR(Y,Phi,z,mun,Sigma,M,A,V)
    Phi <- t(scale(t(Phi)))
    
    ### Update random effects variance
    vn <- nu0 + n
    Dstar <- S0 + t(Phi) %*% (M - A) %*% Phi 
    V <- solve(r2arma::rwishart(vn,solve(Dstar)))
    
    z <- update_z_spot_MCAR(z,Y,Phi,mun,Sigma,pi,1:K)
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
                   SIGMA = SIGMA,
                   K = K,
                   Z = Z,
                   z = z_map)
  return(ret_list)
}