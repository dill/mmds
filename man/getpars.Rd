\name{getpars}
\alias{getpars}
\title{Grab parameter values}
\description{Extract parameter values and create a named list.
}
\usage{
getpars(fpar,mix.terms,zdim=0,z=NULL)
}
\arguments{
   \item{fpar}{The \code{$par} element \code{\link{ds.mixture}} object}
   \item{mix.terms}{Number of mixture components}
   \item{zdim}{Number of covariates.}
   \item{z}{Covariate matrix.}
}
\details{
}
\value{
  result: a named list with elements \code{$key.scale} (giving the key scales) and \code{$mix.prop} giving the mixture proportions. 
}
\references{ 
Miller, D.L. and L. Thomas (in prep.). Mixture model distance sampling detection functions.
 }
\author{David Lawrence Miller}
\seealso{\code{\link{summary.ds.mixture}}, \code{\link{step.ds.mixture}},\code{\link{plot.ds.mixture}}, \code{\link{sim.mix}}, \code{\link{mmds.gof}},\code{\link{ds.mixture}}, \code{\link{mmds.gof}}, \code{\link{fitmix}}, \code{\link{mmds.pars}}, \code{\link{step.ds.mixture}}}
\keyword{Statistical Models, Distance Sampling}
