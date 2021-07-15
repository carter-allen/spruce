// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mvrnormArma
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP _spruce_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// sub_mat
arma::mat sub_mat(arma::vec A, arma::mat X);
RcppExport SEXP _spruce_sub_mat(SEXP ASEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(sub_mat(A, X));
    return rcpp_result_gen;
END_RCPP
}
// col_sum
arma::vec col_sum(arma::mat X);
RcppExport SEXP _spruce_col_sum(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(col_sum(X));
    return rcpp_result_gen;
END_RCPP
}
// update_phi_spot_MCAR
arma::mat update_phi_spot_MCAR(NumericMatrix Y, arma::mat Phi, NumericVector z, List mun, List Sigma, NumericMatrix M, arma::mat A, arma::mat L);
RcppExport SEXP _spruce_update_phi_spot_MCAR(SEXP YSEXP, SEXP PhiSEXP, SEXP zSEXP, SEXP munSEXP, SEXP SigmaSEXP, SEXP MSEXP, SEXP ASEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< List >::type mun(munSEXP);
    Rcpp::traits::input_parameter< List >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(update_phi_spot_MCAR(Y, Phi, z, mun, Sigma, M, A, L));
    return rcpp_result_gen;
END_RCPP
}
// rtn
double rtn(double a, double b, double mu, double sig);
RcppExport SEXP _spruce_rtn(SEXP aSEXP, SEXP bSEXP, SEXP muSEXP, SEXP sigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    rcpp_result_gen = Rcpp::wrap(rtn(a, b, mu, sig));
    return rcpp_result_gen;
END_RCPP
}
// update_t
NumericVector update_t(NumericVector ts, int k, NumericVector k_inds, double Ak, arma::vec xnk, arma::mat Sigmak, NumericMatrix Y, arma::vec munk, double a, double b);
RcppExport SEXP _spruce_update_t(SEXP tsSEXP, SEXP kSEXP, SEXP k_indsSEXP, SEXP AkSEXP, SEXP xnkSEXP, SEXP SigmakSEXP, SEXP YSEXP, SEXP munkSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ts(tsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k_inds(k_indsSEXP);
    Rcpp::traits::input_parameter< double >::type Ak(AkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xnk(xnkSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigmak(SigmakSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type munk(munkSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(update_t(ts, k, k_inds, Ak, xnk, Sigmak, Y, munk, a, b));
    return rcpp_result_gen;
END_RCPP
}
// nnk
int nnk(NumericVector zs, NumericMatrix A, int k, int i);
RcppExport SEXP _spruce_nnk(SEXP zsSEXP, SEXP ASEXP, SEXP kSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(nnk(zs, A, k, i));
    return rcpp_result_gen;
END_RCPP
}
// update_z
NumericVector update_z(NumericVector zs, NumericMatrix Y, List mun, List Sigma, NumericVector pi, NumericVector classes);
RcppExport SEXP _spruce_update_z(SEXP zsSEXP, SEXP YSEXP, SEXP munSEXP, SEXP SigmaSEXP, SEXP piSEXP, SEXP classesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type mun(munSEXP);
    Rcpp::traits::input_parameter< List >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type classes(classesSEXP);
    rcpp_result_gen = Rcpp::wrap(update_z(zs, Y, mun, Sigma, pi, classes));
    return rcpp_result_gen;
END_RCPP
}
// update_z_PG
NumericVector update_z_PG(NumericVector zs, NumericMatrix Y, List mun, List Sigma, NumericMatrix Pi, NumericVector classes);
RcppExport SEXP _spruce_update_z_PG(SEXP zsSEXP, SEXP YSEXP, SEXP munSEXP, SEXP SigmaSEXP, SEXP PiSEXP, SEXP classesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type mun(munSEXP);
    Rcpp::traits::input_parameter< List >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type classes(classesSEXP);
    rcpp_result_gen = Rcpp::wrap(update_z_PG(zs, Y, mun, Sigma, Pi, classes));
    return rcpp_result_gen;
END_RCPP
}
// update_z_MSN
NumericVector update_z_MSN(NumericVector zs, NumericMatrix Y, NumericVector t, List mun, List xin, List Sigma, NumericVector pi, NumericVector classes);
RcppExport SEXP _spruce_update_z_MSN(SEXP zsSEXP, SEXP YSEXP, SEXP tSEXP, SEXP munSEXP, SEXP xinSEXP, SEXP SigmaSEXP, SEXP piSEXP, SEXP classesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< List >::type mun(munSEXP);
    Rcpp::traits::input_parameter< List >::type xin(xinSEXP);
    Rcpp::traits::input_parameter< List >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type classes(classesSEXP);
    rcpp_result_gen = Rcpp::wrap(update_z_MSN(zs, Y, t, mun, xin, Sigma, pi, classes));
    return rcpp_result_gen;
END_RCPP
}
// update_z_spot_MCAR
NumericVector update_z_spot_MCAR(NumericVector zs, NumericMatrix Y, NumericMatrix Phi, List mun, List Sigma, NumericVector pi, NumericVector classes);
RcppExport SEXP _spruce_update_z_spot_MCAR(SEXP zsSEXP, SEXP YSEXP, SEXP PhiSEXP, SEXP munSEXP, SEXP SigmaSEXP, SEXP piSEXP, SEXP classesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< List >::type mun(munSEXP);
    Rcpp::traits::input_parameter< List >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type classes(classesSEXP);
    rcpp_result_gen = Rcpp::wrap(update_z_spot_MCAR(zs, Y, Phi, mun, Sigma, pi, classes));
    return rcpp_result_gen;
END_RCPP
}
// update_z_smooth
NumericVector update_z_smooth(NumericVector zs, NumericMatrix Y, List mun, List Sigma, NumericVector pis, NumericVector classes, double gamma, NumericMatrix M, NumericMatrix A);
RcppExport SEXP _spruce_update_z_smooth(SEXP zsSEXP, SEXP YSEXP, SEXP munSEXP, SEXP SigmaSEXP, SEXP pisSEXP, SEXP classesSEXP, SEXP gammaSEXP, SEXP MSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type mun(munSEXP);
    Rcpp::traits::input_parameter< List >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pis(pisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type classes(classesSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(update_z_smooth(zs, Y, mun, Sigma, pis, classes, gamma, M, A));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spruce_mvrnormArma", (DL_FUNC) &_spruce_mvrnormArma, 3},
    {"_spruce_sub_mat", (DL_FUNC) &_spruce_sub_mat, 2},
    {"_spruce_col_sum", (DL_FUNC) &_spruce_col_sum, 1},
    {"_spruce_update_phi_spot_MCAR", (DL_FUNC) &_spruce_update_phi_spot_MCAR, 8},
    {"_spruce_rtn", (DL_FUNC) &_spruce_rtn, 4},
    {"_spruce_update_t", (DL_FUNC) &_spruce_update_t, 10},
    {"_spruce_nnk", (DL_FUNC) &_spruce_nnk, 4},
    {"_spruce_update_z", (DL_FUNC) &_spruce_update_z, 6},
    {"_spruce_update_z_PG", (DL_FUNC) &_spruce_update_z_PG, 6},
    {"_spruce_update_z_MSN", (DL_FUNC) &_spruce_update_z_MSN, 8},
    {"_spruce_update_z_spot_MCAR", (DL_FUNC) &_spruce_update_z_spot_MCAR, 7},
    {"_spruce_update_z_smooth", (DL_FUNC) &_spruce_update_z_smooth, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_spruce(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
