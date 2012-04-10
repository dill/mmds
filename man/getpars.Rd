\name{getpars}
\alias{getpars}
\title{Grab parameter values}
\usage{
  getpars(fpar, mix.terms, zdim = 0, z = NULL)
}
\arguments{
  \item{fpar}{The \code{$par} element
  \code{\link{ds.mixture}} object.}

  \item{mix.terms}{Number of mixture components.}

  \item{zdim}{Number of covariates.}

  \item{z}{Covariate matrix.}
}
\value{
  a named list with elements \code{$key.scale} (giving the
  key scales) and \code{$mix.prop} giving the mixture
  proportions.
}
\description{
  Extract parameter values and create a named list.
}
\author{
  David L. Miller
}

