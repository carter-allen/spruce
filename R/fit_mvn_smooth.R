#' Spatial multivariate normal mixture model clustering
#'
#' Implement Gibbs sampling for MVN model with spatial smoothing
#'
#' @param Y An n x g matrix of gene expression values. n is the number of cell spots and g is the number of features.
#' @param coords_df An n x 2 data frame or matrix of 2d spot coordinates.  
#' @param K The number of mixture components to fit. 
#' @param r Empirical spatial smoothing
#' @param nsim Number of total MCMC iterations to run.
#' @param burn Number of MCMC iterations to discard as burn in. The number of saved samples is nsim - burn.
#' @param z_init Optional initialized allocation vector. Randomly initialized if NULL. 
#' @param verbose Logical for printing cluster allocations at each iteration.
#'
#' @return a list of posterior samples
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats cov kmeans
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
#' g <- 3 # number of features
#' K <- length(unique(coords_df_sim$z)) # number of clusters (mixture components)
#' pi <- table(z)/length(z) # cluster membership probability
#' 
#' # Cluster Specific Parameters
#' # cluster specific means
#' Mu <- list(
#'   Mu1 = rnorm(g,-2,1),
#'   Mu2 = rnorm(g,-1,1),
#'   Mu3 = rnorm(g,1,1),
#'   Mu4 = rnorm(g,2,1)
#' )
#' # cluster specific variance-covariance
#' S <- matrix(0.5,nrow = g,ncol = g) # y covariance matrix
#' diag(S) <- 1
#' Sig <- list(
#'   Sig1 = S,
#'   Sig2 = S, 
#'   Sig3 = S,
#'   Sig4 = S
#' )
#' 
#' Y <- matrix(0, nrow = n, ncol = g)
#' for(i in 1:n)
#' {
#'   Y[i,] <- mvtnorm::rmvnorm(1,mean = Mu[[z[i]]],sigma = Sig[[z[i]]])
#' }
#' 
#' # sometimes helps to initialize using heuristic like kmeans
#' fitk <- stats::kmeans(Y,4)
#' z_km <- remap_canonical2(fitk$cluster)
#' 
#' # fit model
#' # use more iterations in practice
#' fit1 <- fit_mvn_smooth(Y,coords_df,4,2,10,0,z_km)}

fit_mvn_smooth <- function(Y,
                           coords_df,
                           K,
                           r,
                           nsim = 2000,
                           burn = 1000,
                           z_init = NULL,
                           verbose = FALSE)
{
  # parameters
  n <- nrow(Y) # number of observations
  p <- ncol(Y) # number of features
  pi <- rep(1/K,K) # cluster membership probability
  if(is.null(z_init)) # initialize z
  {
    fit_kmeans <- kmeans(Y,centers = K)
    z_init <- fit_kmeans$cluster
    z <- z_init
  }
  else # user provided initialization
  {
    z <- z_init
    pi <- table(z)/n
  }
  
  # adjacency matrix
  A <- build_knn_graph(coords_df, k = 4)
  m <- colSums(A)
  M <- diag(m)
  
  # priors - shared across clusters
  mu0 <- colMeans(Y)
  L0 <- S0 <- diag(p)
  nu0 <- 2

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
      Sigma[[k]] <- solve(rwishart(nu0+nk, solve(Sn[[k]])))
    }
    
    ### Update cluster indicators
    z <- update_z_smooth(z,Y,mun,Sigma,pi,1:K,r,M,A)
    # remap to address label switching
    z <- remap_canonical2(z)
    pi <- update_props(z,K)
    if(verbose)
    {
      print(update_counts(z,K))
    }
    if(any(update_counts(z,K) < 20))
    {
      z = sample(1:K, size = n, replace = TRUE, prob = rep(1/K,K))
      pi <- update_props(z,K)
    }
    
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
                   W = NULL,
                   coords_df = coords_df,
                   MU = MU,
                   XI = NULL,
                   SIGMA = SIGMA,
                   DELTA = NULL,
                   K = K,
                   Z = Z,
                   z = z_map,
                   z_init = z_init)
  return(ret_list)
}