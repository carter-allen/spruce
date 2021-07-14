// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

// update for traditional multivariate Normal mixture
// [[Rcpp::export]]
NumericVector update_z(NumericVector zs, 
                       NumericMatrix Y,
                       List mun, 
                       List Sigma,
                       NumericVector pi,
                       NumericVector classes) 
{
  int n = zs.length();
  int K = classes.length();
  NumericVector z_ret(n);
  NumericMatrix PZ(2,n);
  for(int i = 0; i < n; i ++)
  {
    NumericVector pj(K);
    for(int k = 0; k < K; k++)
    {
      NumericVector munk = mun[k];
      arma::mat Sigmak = Sigma[k];
      arma::rowvec yj = Y(i,_);
      arma::vec pjk = dmvnrm_arma_fast(yj,munk,Sigmak);
      pj[k] = pjk[0];
    }
    pj = pi*pj / sum(pi*pj);
    //Rcout << i << ":" << pi_star << std::endl;
    NumericVector zi = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pj);
    int z_i = zi[0];
    z_ret[i] = z_i;
  }
  return z_ret;
}

// update for multivariate Normal mixture
// PG multinomial regression w/ no REs
// [[Rcpp::export]]
NumericVector update_z_PG(NumericVector zs, 
                          NumericMatrix Y,
                          List mun, 
                          List Sigma,
                          NumericMatrix Pi,
                          NumericVector classes) 
{
  int n = zs.length();
  int K = classes.length();
  NumericVector z_ret(n);
  NumericMatrix PZ(2,n);
  for(int i = 0; i < n; i ++)
  {
    NumericVector pj(K);
    for(int k = 0; k < K; k++)
    {
      NumericVector munk = mun[k];
      arma::mat Sigmak = Sigma[k];
      arma::rowvec yj = Y(i,_);
      arma::vec pjk = dmvnrm_arma_fast(yj,munk,Sigmak);
      pj[k] = pjk[0];
    }
    NumericVector pii = Pi(i,_);
    NumericVector pjn = pii*pj / sum(pii*pj);
    NumericVector zi = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pjn);
    int z_i = zi[0];
    z_ret[i] = z_i;
  }
  return z_ret;
}

// update for MSN mixture
// [[Rcpp::export]]
NumericVector update_z_MSN(NumericVector zs, 
                           NumericMatrix Y,
                           NumericVector t,
                           List mun, 
                           List xin,
                           List Sigma,
                           NumericVector pi,
                           NumericVector classes) 
{
  int n = zs.length();
  int K = classes.length();
  NumericVector z_ret(n);
  for(int i = 0; i < n; i ++)
  {
    NumericVector pj(K);
    arma::rowvec yj = Y(i,_);
    double ti = t[i];
    for(int k = 0; k < K; k++)
    {
      NumericVector munk = mun[k];
      NumericVector xink = xin[k];
      NumericVector etaik = munk + ti * xink;
      arma::mat Sigmak = Sigma[k];
      arma::vec pjk = dmvnrm_arma_fast(yj,etaik,Sigmak);
      pj[k] = pjk[0];
    }
    pj = pi*pj / sum(pi*pj);
    //Rcout << i << ":" << pi_star << std::endl;
    NumericVector zi = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pj);
    z_ret[i] = zi[0];
  }
  return z_ret;
}