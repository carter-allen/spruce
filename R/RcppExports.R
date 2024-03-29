# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rwishart <- function(df, S) {
    .Call('_spruce_rwishart', PACKAGE = 'spruce', df, S)
}

mvrnormArma <- function(n, mu, sigma) {
    .Call('_spruce_mvrnormArma', PACKAGE = 'spruce', n, mu, sigma)
}

sub_mat <- function(A, X) {
    .Call('_spruce_sub_mat', PACKAGE = 'spruce', A, X)
}

col_sum <- function(X) {
    .Call('_spruce_col_sum', PACKAGE = 'spruce', X)
}

update_phi_spot_MCAR <- function(Y, Phi, z, mun, Sigma, M, A, L) {
    .Call('_spruce_update_phi_spot_MCAR', PACKAGE = 'spruce', Y, Phi, z, mun, Sigma, M, A, L)
}

rtn <- function(a, b, mu, sig) {
    .Call('_spruce_rtn', PACKAGE = 'spruce', a, b, mu, sig)
}

update_t <- function(ts, k, k_inds, Ak, xnk, Sigmak, Y, munk, a, b) {
    .Call('_spruce_update_t', PACKAGE = 'spruce', ts, k, k_inds, Ak, xnk, Sigmak, Y, munk, a, b)
}

nnk <- function(zs, A, k, i) {
    .Call('_spruce_nnk', PACKAGE = 'spruce', zs, A, k, i)
}

update_counts <- function(zs, K0) {
    .Call('_spruce_update_counts', PACKAGE = 'spruce', zs, K0)
}

update_props <- function(zs, K0) {
    .Call('_spruce_update_props', PACKAGE = 'spruce', zs, K0)
}

update_z <- function(zs, Y, mun, Sigma, pi, classes) {
    .Call('_spruce_update_z', PACKAGE = 'spruce', zs, Y, mun, Sigma, pi, classes)
}

update_z_PG <- function(zs, Y, mun, Sigma, Pi, classes) {
    .Call('_spruce_update_z_PG', PACKAGE = 'spruce', zs, Y, mun, Sigma, Pi, classes)
}

update_z_PG_smooth <- function(zs, Y, mun, Sigma, Pi, classes, gamma, M, A) {
    .Call('_spruce_update_z_PG_smooth', PACKAGE = 'spruce', zs, Y, mun, Sigma, Pi, classes, gamma, M, A)
}

update_z_MSN_PG_smooth <- function(zs, Y, t, mun, xin, Sigma, Pi, classes, gamma, M, A) {
    .Call('_spruce_update_z_MSN_PG_smooth', PACKAGE = 'spruce', zs, Y, t, mun, xin, Sigma, Pi, classes, gamma, M, A)
}

update_z_MSN <- function(zs, Y, t, mun, xin, Sigma, pi, classes) {
    .Call('_spruce_update_z_MSN', PACKAGE = 'spruce', zs, Y, t, mun, xin, Sigma, pi, classes)
}

update_z_MSN_smooth <- function(zs, Y, t, mun, xin, Sigma, pi, classes, gamma, M, A) {
    .Call('_spruce_update_z_MSN_smooth', PACKAGE = 'spruce', zs, Y, t, mun, xin, Sigma, pi, classes, gamma, M, A)
}

update_z_spot_MCAR <- function(zs, Y, Phi, mun, Sigma, pi, classes) {
    .Call('_spruce_update_z_spot_MCAR', PACKAGE = 'spruce', zs, Y, Phi, mun, Sigma, pi, classes)
}

update_z_spot_PG_MCAR <- function(zs, Y, Phi, mun, Sigma, Pi, classes) {
    .Call('_spruce_update_z_spot_PG_MCAR', PACKAGE = 'spruce', zs, Y, Phi, mun, Sigma, Pi, classes)
}

update_z_spot_PG_MCAR_smooth <- function(zs, Y, Phi, mun, Sigma, Pi, classes, gamma, M, A) {
    .Call('_spruce_update_z_spot_PG_MCAR_smooth', PACKAGE = 'spruce', zs, Y, Phi, mun, Sigma, Pi, classes, gamma, M, A)
}

update_z_smooth <- function(zs, Y, mun, Sigma, pis, classes, gamma, M, A) {
    .Call('_spruce_update_z_smooth', PACKAGE = 'spruce', zs, Y, mun, Sigma, pis, classes, gamma, M, A)
}

