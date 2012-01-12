\name{reparam.pi}
\alias{inv.reparam.pi}
\alias{reparam.pi}
\title{Reparametrize mixture proportions}
\description{Re-parameterise the mixture proportions so that when there is more than a 2-point mixture, the proportions sum to 1.
}
\usage{
inv.reparam.pi(mix.prop,lastpar=FALSE)
reparam.pi(thetas)
}
\arguments{
   \item{mix.prop}{A set of mixture proportions}
   \item{lastpar}{Is the last parameter provided, i.e. does sum(mix.prop)=1?}
   \item{thetas}{Mixture proportions in their parametrisation for optimization}
}
\details{See Miller and Thomas for information on exactly how these are calculated. Thanks go to David Borchers for proposing the trick.
}
\value{
  result: a vector of parameters in one or other parametrisation 
}
\references{ 
Miller, D.L. and L. Thomas (in prep.). Mixture model distance sampling detection functions.
 }
\author{David Lawrence Miller}
\seealso{\code{\link{summary.ds.mixture}}, \code{\link{step.ds.mixture}},\code{\link{plot.ds.mixture}}, \code{\link{sim.mix}}, \code{\link{mmds.gof}},\code{\link{ds.mixture}}, \code{\link{mmds.gof}}, \code{\link{fitmix}}, \code{\link{mmds.pars}}, \code{\link{step.ds.mixture}}}
\examples{
library(mmds)
reparam.pi(inv.reparam.pi(0.3))
reparam.pi(inv.reparam.pi(c(0.3,0.4,0.1),lastpar=TRUE))
}
\keyword{Statistical Models, Distance Sampling}
