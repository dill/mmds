\name{ds.mixture}
\alias{ds.mixture}
\title{A fitted Mixture Model Detection Function Object}
\description{
  The fitted mixture model detection function object
  returned by \code{\link{fitmix}}. Knowledge of most of
  this is not useful. Use \code{link{summary.ds.mixture}}
  for result summaries.
}
\section{Details}{
  A \code{\link{ds.mixture}} object has the following
  elements: \tabular{ll}{ distance \tab Vector of distances
  used in the analysis.\cr likelihood \tab Value of the
  log-likelihood at the maxima.\cr pars \tab Parmeter
  estimates. See \code{\link{mmds.pars}} for more
  information.\cr mix.terms \tab Number of mixture terms
  fit.\cr width \tab Truncation distance used.\cr z \tab
  List containing the matrix of covariates used. Output
  from \code{\link{model.matrix}}.\cr zdim \tab Number of
  columns of \code{z}. See \code{\link{mmds.pars}} for more
  information.\cr hessian \tab Hessian matrix at the
  maxima.\cr pt \tab Logical indicating whether the data
  were from a point transect survey.\cr data \tab Data
  frame after truncation.\cr ftype \tab Type of detection
  function.\cr ctrl.options \tab Options passed to
  \code{\link{optim}}.\cr showit \tab Debug level.\cr
  opt.method \tab Optimisation method used.\cr usegrad \tab
  Were analytic gradients used?\cr model.formula \tab Model
  formula.\cr mu \tab Per-observation effective trip
  width/effective area of detection.\cr pa.vec \tab Vector
  of per-observation detectabilities.\cr N \tab Estimate of
  N in the covered area (Horvitz-Thompson).\cr pa \tab
  Average detectability.\cr pars.se \tab Standard errors of
  the parameters.\cr N.se \tab Standard error of the
  Horvitz-Thompson estimate of the abundance.\cr pa.se \tab
  Standard error of the average detectability.\cr aic \tab
  AIC of the fitted model.\cr cvm \tab Cramer-von Mises GoF
  test results. List containing: \code{p}, the p-value and
  \code{W}, the test statistic.\cr ks \tab
  Kolmogorov-Smirnov test results. List containing:
  \code{p}, the p-value and \code{Dn}, the test statistic.
  See \code{\link{mmds.gof}} for more information.\cr}
}

\section{Note}{
  \code{ds.mixture} objects can be passed to
  \code{\link{step.ds.mixture}} to select number of mixture
  components based on AIC score.
}
\author{
  David L. Miller
}

