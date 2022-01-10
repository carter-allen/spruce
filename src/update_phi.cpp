#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// adapted from r2arma package (https://github.com/coatless-rpkg/r-to-armadillo/blob/master/src/distributions.cpp)
// [[Rcpp::export]]
arma::mat rwishart(unsigned int df, const arma::mat& S)
{
  // Dimension of returned wishart
  unsigned int m = S.n_rows;
  
  // Z composition:
  // sqrt chisqs on diagonal
  // random normals below diagonal
  // misc above diagonal
  arma::mat Z(m,m);
  
  // Fill the diagonal
  for(unsigned int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i));
  }
  
  // Fill the lower matrix with random guesses
  for(unsigned int j = 0; j < m; j++){  
    for(unsigned int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1);
    }
  }
  
  // Lower triangle * chol decomp
  arma::mat C = arma::trimatl(Z).t() * arma::chol(S);
  
  // Return random wishart
  return C.t()*C;
}


// Rcpp for sampling multivariate normal
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::mat sub_mat(arma::vec A, arma::mat X)
{
  int P = X.n_cols;
  int n = X.n_rows;
  int n_a = arma::accu(A);
  arma::mat X_ret(n_a,P);
  int n_i = 0;
  for(int i = 0; i < n; i++)
  {
    if(A[i] == 1)
    {
      X_ret.row(n_i) = X.row(i);
      n_i = n_i + 1;
    }
  }
  return(X_ret);
}

// [[Rcpp::export]]
arma::vec col_sum(arma::mat X)
{
  int P = X.n_cols;
  arma::vec s_ret(P);
  for(int p = 0; p < P; p++)
  {
    s_ret[p] = arma::accu(X.col(p));
  }
  return(s_ret);
}

// Rcpp for updating MCAR random intercepts
// Non-cluster specific
// Phi is n x p
// [[Rcpp::export]]
arma::mat update_phi_spot_MCAR(NumericMatrix Y,
                               arma::mat Phi,
                               NumericVector z, 
                               List mun, 
                               List Sigma,
                               NumericMatrix M,
                               arma::mat A,
                               arma::mat L) {
  int n = z.length();
  int p = Y.ncol();
  arma::mat Phi_ret = Phi;
  arma::mat Vi(p,p);
  arma::vec mi(p);
  arma::vec mui(p);
  for(int i = 0; i < n; i ++)
  {
    arma::vec yi = Y(i,_);
    arma::vec Ai = A.col(i);
    int mi = M(i,i);
    int zi = z[i];
    arma::mat Sigmak = Sigma[zi - 1];
    arma::vec muzi = mun[zi - 1];
    Vi = inv(inv(Sigmak) + mi * inv(L));
    mui = Vi * (inv(Sigmak) * (yi - muzi) + inv(L) * col_sum(sub_mat(Ai,Phi)));
    Phi_ret.row(i) = mvrnormArma(1,mui,Vi);
    Phi_ret.row(i) = Phi_ret.row(i) - mean(Phi_ret.row(i));
  }
  return Phi_ret;
}
