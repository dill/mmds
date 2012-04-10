\name{mmds.pars}
\alias{mmds.pars}
\title{Parameters in \code{mmds}}
\description{
  The internal parametrisation used in \code{mmds} is not
  directly interpretable. This man page aims to explain how
  to interpret the parameters and transform them into
  useful information.
}
\section{Details}{
  Parametrisation works differently for the scale
  parameters of the half-normals and for the mixture
  proportions.

  The scale parameters of the half-normals (or their
  constituent parameters in the case of a covariate model)
  are given on the log scale.

  The mixture proportions are transformed to a
  parametrisation that allows values to lie over all of the
  real line (see Miller and Thomas for details).

  The parameter vector is made up of the scale parameters
  followed by the mixture parameters. In the non-covariate
  case the former is the length of the number of mixtures
  (\code{mix.terms}) and the latter is of length
  \code{mix.terms}-1. When the model has covariates the
  scale paramters are given as \code{mix.terms} intercepts
  followed by the covariate parameters.

  The function \code{\link{getpars}} transforms the
  parameters (\code{$pars} element) in a
  \code{\link{ds.mixture}} object to a named list.

  Calling \code{\link{summary.ds.mixture}} will show the
  mixture proportions.

  The parameter \code{initialvalues} supplied to
  \code{\link{fitmix}}.
}
\author{
  David L. Miller
}

