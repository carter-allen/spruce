#' Multivariate normal spatial mixture model clustering w/ PG multinomial regression on membership probabilities with spatial smoothing
#'
#' Implement Gibbs sampling for MVN model with MCAR spatial random effects w/ PG multinomial regression on membership probabilities and CAR random ints in multinomial regression model with spatial smoothing.
#'
#' @param Y An n x g matrix of gene expression values. n is the number of cell spots and g is the number of features.
#' @param W An n x v matrix of covariates to predict cluster membership. Should include an intercept (i.e., first column is 1)
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
#' @importFrom stats cov rnorm rgamma kmeans
#' @importFrom MCMCpack rdirichlet

fit_mvn_PG_CAR_MCAR_smooth <- function(Y,
                                       W,
                                       coords_df,
                                       K,
                                       r = 3,
                                       nsim = 2000,
                                       burn = 1000,
                                       z_init = NULL, 
                                       verbose = FALSE)
{
  # parameters
  n <- nrow(coords_df) # number of observations
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
  
  # adjacency matrix
  A <- build_knn_graph(coords_df, k = 4)
  m <- colSums(A)
  M <- diag(m)
  Q <- M - A
  
  # random effects
  PHI <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  PSI <- matrix(0, nrow = nrow(Y), ncol = K-1)
  
  # priors - shared across clusters
  mu0 <- colMeans(Y)
  L0 <- S0 <- V <- diag(p)
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
  Nu2 <- rep(1,K-1)
  
  # Empty sample storage
  MU <- SIGMA <- vector("list",K)
  
  n_save <- nsim - burn
  Z <- matrix(0,nrow = n_save,ncol = n)
  Phi <- PHI
  
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
      Ybar[[k]] <- colMeans(Y[z == k,] - Phi[z == k,])
      
      ### update mu - cluster specific
      Ln[[k]] <- solve(solve(L0) + nk*solve(Sigma[[k]]))
      mn[[k]] <- Ln[[k]] %*% (solve(L0) %*% mu0 + nk*solve(Sigma[[k]]) %*% Ybar[[k]]) 
      mun[[k]] <- mvrnormArma(1 ,mn[[k]],Ln[[k]])
      
      ### update Sigma - cluster specific 
      Sn[[k]] <- S0 + (t(Y[z == k,]) - t(Phi[z == k,]) - c(mun[[k]])) %*% t(t(Y[z == k,]) - t(Phi[z == k,]) - c(mun[[k]])) 
      Sigma[[k]] <- solve(rwishart(nu0+nk, solve(Sn[[k]])))
    }
    
    ### Update random effects Phi
    Phi <- update_phi_spot_MCAR(Y,Phi,z,mun,Sigma,M,A,V)
    Phi <- t(scale(t(Phi)))
    
    ### Update random effects variance
    vn <- nu0 + n
    Dstar <- S0 + t(Phi) %*% (M - A) %*% Phi 
    V <- solve(rwishart(vn,solve(Dstar)))
    
    z <- update_z_spot_PG_MCAR_smooth(z,Y,Phi,mun,Sigma,PI,1:K,r,M,A)
    # remap to address label switching
    z <- remap_canonical2(z)
    
    # Update multinomial regression parameters
    W <- as.matrix(W) # enforce W is a matrix
    for(k in 1:(K-1))
    {
      deltak <- delta[,k]
      deltanotk <- delta[,-k]
      PSIk <- PSI[,k]
      PSInotk <- PSI[,-k]
      uk <- 1*(z == (k+1))
      ck <- log(1 + rowSums(exp(W %*% deltanotk + PSInotk)))
      eta <- W %*% deltak + PSIk - ck
      w <- rpg(n, 1, eta)
      ukstr <- (uk - 1/2)/w + ck
      
      D <- solve(D0 + crossprod(W*sqrt(w)))  
      d <- D %*% (D0 %*% delta0 + t(w*W) %*% ukstr)
      deltak <- c(rmvnorm(1,d,D))
      delta[,k] <- deltak
      
      Mk <- (get_psi_sums(PSIk,A)/m + ukstr - W%*%deltak)/(m^2/Nu2[k] + 1/w^2)
      Vk <- 1/(m^2/Nu2[k] + 1/w^2)
      PSI[,k] <- rnorm(n = n, mean = Mk, sd = sqrt(Vk))
      
      Nu2[k] <- 1/rgamma(1,.01 + (n-1)/2,.01 + (t(PSI[,k]) %*% Q %*% PSI[,k])/2)
    }
    # delta <- Delta
    eta <- cbind(rep(0,n),W%*%delta + PSI)
    PI <- exp(eta)/(1+apply(as.matrix(exp(eta[,-1])),1,sum))
    pi <- table(z)/n
    
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
                   coords_df = coords_df,
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