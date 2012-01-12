\name{mmds.pars}
\alias{mmds.pars}
\title{Parameters in \code{mmds}}
\description{
The internal parametrisation used in \code{mmds} is not directly interpretable. This man page aims to explain how to interpret the parameters and transform them into useful information.
}
\details{
Parametrisation works differently for the scale parameters of the half-normals and for the mixture proportions.

The scale parameters of the half-normals (or their constituent parameters in the case of a covariate model) are given on the log scale.

The mixture proportions are transformed to a parametrisation that allows values to lie over all of the real line (see Miller and Thomas for details).

The parameter vector is made up of the scale parameters followed by the mixture parameters. In the non-covariate case the former is the length of the number of mixtures (\code{mix.terms}) and the latter is of length \code{mix.terms}-1. When the model has covariates the scale paramters are given as \code{mix.terms} intercepts followed by the covariate parameters.

The function \code{\link{getpars}} transforms the parameters (\code{$pars} element) in a \code{\link{ds.mixture}} object to a named list.

Calling \code{\link{summary.ds.mixture}} will show the mixture proportions.

The parameter \code{initialvalues} supplied to \code{\link{fitmix}}.
}
\references{ 
Miller, D.L. and L. Thomas (in prep.). Mixture model distance sampling detection functions.
 }
\author{David Lawrence Miller}
\note{
\code{ds.mixture} objects can be passed to \code{\link{step.ds.mixture}} to select number of mixture components based on AIC score. 
}
\seealso{\code{\link{summary.ds.mixture}}, \code{\link{step.ds.mixture}},\code{\link{plot.ds.mixture}}, \code{\link{sim.mix}}, \code{\link{mmds.gof}},\code{\link{ds.mixture}}, \code{\link{mmds.gof}}, \code{\link{fitmix}}, \code{\link{mmds.pars}}, \code{\link{step.ds.mixture}}}
\keyword{Statistical Models, Distance Sampling}
