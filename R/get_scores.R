#' Calculate cluster uncertainty
#'
#' Use posterior estimates to calculate uncertainty scores
#'
#' @param fit A model fit returned by one of the fit_*_PG model functions
#'
#' @return An n x (K + 1) matrix. First K columns are continuous phenotypes, and last column is uncertainty scores 
#' @export
#' @importFrom mvtnorm dmvnorm
#' @examples 
#' # scores_df <- get_scores(fit)
get_scores <- function(fit)
{
  Y <- fit$Y
  n <- nrow(Y)
  W <- fit$W
  v <- ncol(W)
  MU <- fit$MU
  K <- fit$K
  p <- ncol(Y)
  SIGMA <- fit$SIGMA
  DELTA <- fit$DELTA
  z_map <- fit$z
  
  delta_post <- matrix(colMeans(DELTA),nrow = v,ncol = K-1)
  eta_post <- cbind(rep(0,n),W%*%delta_post)
  PI_post <- exp(eta_post)/(1+apply(as.matrix(exp(eta_post[,-1])),1,sum))
  pi_post <- table(z_map)/n
  PI_d_post <- matrix(0,nrow = n,ncol = K)
  
  mu_post <- lapply(MU, colMeans)
  sig_post <- lapply(SIGMA, colMeans)
  
  ro <- matrix(0,nrow = n,ncol = K)
  
  for(i in 1:n)
  {
    pj <- rep(0,K)
    yi <- Y[i,]
    for(k in 1:K)
    {
      mu_k_post <- mu_post[[k]]
      sigma_k_post <- matrix(sig_post[[k]],
                             nrow = p,
                             ncol = p)
      pj[k] <- dmvnorm(yi,mu_k_post,sigma_k_post)
    }
    pii <- PI_post[i,]
    pjn <- pii*pj/sum(pii*pj)
    PI_d_post[i,] <- pjn
  }
  
  u_score <- rep(0,n)
  for(i in 1:n)
  {
    u_score[i] <- 1-max(PI_d_post[i,])
  }
  
  ret <- cbind(PI_d_post,u_score)
  colnames(ret) <- c(paste0("k",1:K,"_score"),"u_score")
  
  return(ret)
}