% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mvn_PG_CAR_smooth.R
\name{fit_mvn_PG_CAR_smooth}
\alias{fit_mvn_PG_CAR_smooth}
\title{Multivariate normal mixture model clustering - PG multinom regression w/ CAR random effect and spatial smoothing}
\usage{
fit_mvn_PG_CAR_smooth(
  Y,
  W,
  coords_df,
  K,
  r = 3,
  nsim = 2000,
  burn = 1000,
  z_init = NULL
)
}
\arguments{
\item{Y}{An n x g matrix of gene expression values. n is the number of cell spots and g is the number of features.}

\item{W}{An n x v matrix of covariates to predict cluster membership. Should include an intercept (i.e., first column is 1)}

\item{coords_df}{An n x 2 data frame or matrix of 2d spot coordinates.}

\item{K}{The number of mixture components to fit.}

\item{r}{Empirical spatial smoothing}

\item{nsim}{Number of total MCMC iterations to run.}

\item{burn}{Number of MCMC iterations to discard as burn in. The number of saved samples is nsim - burn.}

\item{z_init}{Optional initialized allocation vector. Randomly initialized if NULL.}
}
\value{
a list of posterior samples
}
\description{
Implement Gibbs sampling for MVN model. Includes fixed effects multinomial regression w/ CAR random intercepts on cluster indicators using Polya-Gamma data augmentation and spatial smoothing.
}
