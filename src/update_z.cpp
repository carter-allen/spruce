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

// [[Rcpp::export]]
int nnk(NumericVector zs,
        NumericMatrix A,
        int k,
        int i)
{
  // get neighbors of i
  int n = A.ncol();
  int count = 0;
  for(int j = 0; j < n; j++)
  {
    int zj = zs[j];
    if(A(i,j) == 1 && zj == k)
    {
      count = count + 1;
    }
  }
  return count;
}

// [[Rcpp::export]]
NumericVector update_counts(NumericVector zs, int K0) {
  NumericVector counts(K0);
  int n = zs.length();
  for(int k = 1; k <= K0; k++)
  {
    for(int i = 0; i < n; i++)
    {
      if(zs[i] == k)
      {
        counts[k-1] = counts[k-1] + 1;
      }
    }
  }
  return counts;
}

// [[Rcpp::export]]
NumericVector update_props(NumericVector zs, int K0) {
  NumericVector counts(K0);
  int n = zs.length();
  for(int k = 1; k <= K0; k++)
  {
    for(int i = 0; i < n; i++)
    {
      if(zs[i] == k)
      {
        counts[k-1] = counts[k-1] + 1;
      }
    }
  }
  return counts/n;
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
      arma::vec pjk = dmvnrm_arma_fast(yj,munk,Sigmak,true);
      pj[k] = pjk[0];
    }
    //Rcout << i << ":" << pj << std::endl;
    NumericVector pj2 (pj.length());
    if(any(exp(pj) == 0).is_true())
    {
      pj2[which_max(pj)] = 1;
    }
    else
    {
      pj = exp(pj);
      pj2 = pi*pj / sum(pi*pj);
    }
    //Rcout << i << ":" << pj2 << std::endl;
    NumericVector zi = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pj2);
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

// update for multivariate Normal mixture
// PG multinomial regression w/ spatial smoothing
// [[Rcpp::export]]
NumericVector update_z_PG_smooth(NumericVector zs, 
                                 NumericMatrix Y,
                                 List mun, 
                                 List Sigma,
                                 NumericMatrix Pi,
                                 NumericVector classes,
                                 double gamma,
                                 NumericMatrix M,
                                 NumericMatrix A) 
{
  int n = zs.length();
  int K = classes.length();
  NumericVector z_ret(n);
  for(int i = 0; i < n; i ++)
  {
    NumericVector pj(K);
    arma::rowvec yj = Y(i,_);
    int mi = M(i,i);
    NumericVector nnk_count(K);
    for(int k = 0; k < K; k++)
    {
      NumericVector munk = mun[k];
      arma::mat Sigmak = Sigma[k];
      arma::vec pjk = dmvnrm_arma_fast(yj,munk,Sigmak,true);
      nnk_count[k] = nnk(zs,A,k+1,i);
      pj[k] = pjk[0];
    }
    NumericVector pj2 (pj.length());
    if(any(exp(pj) == 0).is_true())
    {
      pj2[which_max(pj)] = 1;
    }
    else
    {
      pj = exp(pj);
      pj = pj*exp((gamma/mi)*2*nnk_count);
      NumericVector pii = Pi(i,_);
      pj2 = pii*pj / sum(pii*pj);
    }
    //Rcout << i << ":" << pj2 << std::endl;
    NumericVector zi = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pj2);
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

// update for multivariate Normal mixture
// with spot-level MCAR spatial random effects
// non-cluster specific RIs
// [[Rcpp::export]]
NumericVector update_z_spot_MCAR(NumericVector zs, 
                                 NumericMatrix Y,
                                 NumericMatrix Phi,
                                 List mun, 
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
    NumericVector phii = Phi(i,_);
    for(int k = 0; k < K; k++)
    {
      NumericVector munk = mun[k];
      NumericVector etaik = phii + munk;
      arma::mat Sigmak = Sigma[k];
      arma::rowvec yj = Y(i,_);
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

// update for multivariate Normal mixture
// with spot-level MCAR spatial random effects
// non-cluster specific RIs
// multinomial PG regressions
// [[Rcpp::export]]
NumericVector update_z_spot_PG_MCAR(NumericVector zs, 
                                    NumericMatrix Y,
                                    NumericMatrix Phi,
                                    List mun, 
                                    List Sigma,
                                    NumericMatrix Pi,
                                    NumericVector classes) 
{
  int n = zs.length();
  int K = classes.length();
  NumericVector z_ret(n);
  for(int i = 0; i < n; i ++)
  {
    NumericVector pj(K);
    NumericVector phii = Phi(i,_);
    for(int k = 0; k < K; k++)
    {
      NumericVector munk = mun[k];
      NumericVector etaik = phii + munk;
      arma::mat Sigmak = Sigma[k];
      arma::rowvec yj = Y(i,_);
      arma::vec pjk = dmvnrm_arma_fast(yj,etaik,Sigmak);
      pj[k] = pjk[0];
    }
    NumericVector pii = Pi(i,_);
    pj = pii*pj / sum(pii*pj);
    //Rcout << i << ":" << pi_star << std::endl;
    NumericVector zi = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pj);
    z_ret[i] = zi[0];
  }
  return z_ret;
}

// update for MVN mixture
// [[Rcpp::export]]
NumericVector update_z_smooth(NumericVector zs, 
                              NumericMatrix Y,
                              List mun, 
                              List Sigma,
                              NumericVector pis,
                              NumericVector classes,
                              double gamma,
                              NumericMatrix M,
                              NumericMatrix A) 
{
  int n = zs.length();
  int K = classes.length();
  NumericVector z_ret(n);
  for(int i = 0; i < n; i ++)
  {
    NumericVector pj(K);
    arma::rowvec yj = Y(i,_);
    int mi = M(i,i);
    for(int k = 0; k < K; k++)
    {
      NumericVector munk = mun[k];
      arma::mat Sigmak = Sigma[k];
      arma::vec pjk = dmvnrm_arma_fast(yj,munk,Sigmak);
      int nnk_count = nnk(zs,A,k+1,i);
      pj[k] = pjk[0]*exp((gamma/mi)*2*nnk_count);
    }
    pj = pis*pj / sum(pis*pj);
    NumericVector zi = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pj);
    z_ret[i] = zi[0];
  }
  return z_ret;
}