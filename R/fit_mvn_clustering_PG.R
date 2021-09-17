#' Multivariate normal mixture model clustering - PG multinom regression
#'
#' Implement Gibbs sampling for MVN model with no spatial random effects. Includes fixed effects multinomial regression on cluster indicators using Polya-Gamma data augmentation.
#'
#' @param Y An n x g matrix of gene expression values. n is the number of cell spots and g is the number of features.
#' @param W An n x v matrix of covariates to predict cluster membership. Should include an intercept (i.e., first column is 1)
#' @param K The number of mixture components to fit. 
#' @param nsim Number of total MCMC iterations to run.
#' @param burn Number of MCMC iterations to discard as burn in. The number of saved samples is nsim - burn.
#' @param z_init Optional initialized allocation vector. Randomly initialized if NULL. 
#'
#' @return a list of posterior samples
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom mvtnorm rmvnorm
#' @importFrom BayesLogit rpg
#' @importFrom stats cov
#' @examples 
#' \dontrun{
#' # parameters
#' n <- 1000 # number of observations
#' g <- 3 # number of features
#' K <- 3 # number of clusters (mixture components)
#' pi <- rep(1/K,K) # cluster membership probability
#' z <- sample(1:K, size = n, replace = TRUE, prob = pi) # cluster indicators
#' z <- remap_canonical2(z)
#' 
#' # multinomial clustering model
#' w1 <- rbinom(n,1,0.5) # binary predictor
#' W <- cbind(1,w1) # design matrix
#' v <- ncol(W) # number of multinomial predictors
#' Delta <- matrix(truncnorm::rtruncnorm(n = v*(K-1),
#'                                       a = -1, 
#'                                       b = 1,
#'                                       mean = 0),
#'                 nrow = v,
#'                 ncol = K-1) # true multinomial regression coefficients
#' Eta <- W %*% Delta # true eta term for each observation
#' U <- matrix(0,nrow = n,ncol = K) # make empty U matrix
#' U[,1] <- QRM::rGumbel(n, mu = 0, sigma = 1)
#' for(i in 2:K) # populate U matrix with gumbel values
#' {
#'   U[,i] <- QRM::rGumbel(n, mu = Eta[,i-1], sigma = 1)
#' }
#' c_true <- c(apply(U,1,which.max))
#' pi_true <- table(c_true)/n
#' z <- c_true
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
#' fit1 <- fit_mvn_clustering_PG(Y,W,3,10,0)}

fit_mvn_clustering_PG <- function(Y,W,K,nsim = 2000,burn = 1000,z_init = NULL)
{
  # parameters
  n <- nrow(Y) # number of observations
  p <- ncol(Y) # number of features
  v <- ncol(W) # number of multinomial predictors
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
    
    ### Update cluster indicators
    z <- update_z_PG(z,Y,mun,Sigma,PI,1:K)
    # remap to address label switching
    z <- remap_canonical2(z)
    
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
      deltak <- c(rmvnorm(1,d,D))
      delta[,k] <- deltak
    }
    # delta <- Delta
    eta <- cbind(rep(0,n),W%*%delta)
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
                   MU = MU,
                   DELTA = DELTA,
                   SIGMA = SIGMA,
                   K = K,
                   Z = Z,
                   z = z_map)
  return(ret_list)
}