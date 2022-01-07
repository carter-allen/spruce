#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Rcpp for sampling univariate truncated normal
// [[Rcpp::export]]
double rtn (double a, double b, double mu, double sig) {
  double draw = 0.0;
  bool valid = false;
  while (!valid) {
    NumericVector cand1 = rnorm(1, mu, sig);
    if ((cand1(0) <= b) & (cand1(0) >= a)) {
      draw = cand1(0);
      valid = true;}
  }
  return(draw);
}

// Rcpp for updating MSN random intercepts
// [[Rcpp::export]]
NumericVector update_t(NumericVector ts,
                       int k,
                       NumericVector k_inds,
                       double Ak,
                       arma::vec xnk,
                       arma::mat Sigmak,
                       NumericMatrix Y,
                       arma::vec munk,
                       double a, 
                       double b)
{
  int nk = k_inds.length();
  for(int i = 0; i < nk; i++)
  {
    int ik = k_inds[i];
    //Rcout << ik << std::endl;
    arma::vec yik = Y(ik - 1,_);
    double aik = Ak * (xnk.t() * inv(Sigmak) * (yik - munk)).eval()(0,0);
    double tik = rtn(a,b,aik,sqrt(Ak));
    ts[ik - 1] = tik;
  }
  return ts;
}