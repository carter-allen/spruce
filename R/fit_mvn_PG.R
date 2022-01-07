#' Multivariate normal mixture model clustering - PG multinom regression
#'
#' Implement Gibbs sampling for MVN model. Includes fixed effects multinomial regression on cluster indicators using Polya-Gamma data augmentation.
#'
#' @param Y An n x g matrix of gene expression values. n is the number of cell spots and g is the number of features.
#' @param W An n x v matrix of covariates to predict cluster membership. Should include an intercept (i.e., first column is 1)
#' @param K The number of mixture components to fit. 
#' @param nsim Number of total MCMC iterations to run.
#' @param burn Number of MCMC iterations to discard as burn in. The number of saved samples is nsim - burn.
#' @param z_init Optional initialized allocation vector. Initialized with hierarchical clustering if NULL. 
#' @param verbose Logical for printing cluster allocations at each iteration.
#'
#' @return a list of posterior samples
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom mvtnorm rmvnorm
#' @importFrom BayesLogit rpg
#' @importFrom stats cov kmeans

fit_mvn_PG <- function(Y,
                       W,
                       K,
                       nsim = 2000,
                       burn = 1000,
                       z_init = NULL, 
                       verbose = FALSE)
{
  # parameters
  n <- nrow(Y) # number of observations
  p <- ncol(Y) # number of features
  v <- ncol(W) # number of multinomial predictors
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
  
  # priors - shared across clusters
  mu0 <- colMeans(Y)
  L0 <- S0 <- diag(p)
  nu0 <- 2
  delta0 <- rep(0,v) # prior mean for delta coefficients (multinomial regression)
  D0 <- diag(1,v) # prior covariance for delta coefficients (multinomial regression)
  
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
  DELTA <- matrix(0,nrow = n_save,ncol = v*(K-1))
  delta <- matrix(0,nrow = v,ncol = K-1)
  eta <- cbind(rep(0,n),W%*%delta)
  PI <- exp(eta)/(1+apply(as.matrix(exp(eta[,-1])),1,sum))
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
    
    # Update multinomial regression parameters
    W <- as.matrix(W) # enforce W is a matrix
    for(k in 1:(K-1))
    {
      deltak <- delta[,k]
      deltanotk <- delta[,-k]
      uk <- 1*(z == (k+1))
      ck <- log(1 + rowSums(exp(W %*% deltanotk)))
      eta <- W %*% deltak - ck
      w <- rpg(n, 1, eta)
      ukstr <- (uk - 1/2)/w + ck
      
      D <- solve(D0 + crossprod(W*sqrt(w)))  
      d <- D %*% (D0 %*% delta0 + t(w*W) %*% ukstr)
      deltak <- c(mvrnormArma(1,d,D))
      delta[,k] <- deltak
    }
    # delta <- Delta
    eta <- cbind(rep(0,n),W%*%delta)
    PI <- exp(eta)/(1+apply(as.matrix(exp(eta[,-1])),1,sum))
    pi <- update_props(z,K)
    
    ### Update cluster indicators
    z = update_z_PG(z,Y,mun,Sigma,PI,1:K)
    if(verbose)
    {
      print(update_counts(z,K))
    }
    if(any(update_counts(z,K) < 20))
    {
      z = z_init
      pi = update_props(z,K)
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
      DELTA[iter,] <- c(delta)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  run.time<-proc.time()-start.time
  print(paste("Finished MCMC after",run.time[1],"seconds"))
  
  z_map <- apply(Z, 2, get_map)
  
  ret_list <- list(Y = Y,
                   W = W,
                   coords_df = NULL,
                   MU = MU,
                   XI = NULL,
                   SIGMA = SIGMA,
                   DELTA = DELTA,
                   K = K,
                   Z = Z,
                   z = z_map,
                   z_init = z_init)
  return(ret_list)
}