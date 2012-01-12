\name{step.ds.mixture}
\alias{step.ds.mixture}
\title{Stepwise selection of mixture components}
\description{
Uses AIC to select the number of mixture components.
}
\usage{
step.ds.mixture(ds.object,max.terms=4)
}
\arguments{
   \item{ds.object}{\code{\link{ds.mixture}} object.}
   \item{max.terms}{Maximum number of mixture components to attempt to fit. Default 4.}
}
\details{
This routine is most useful during model building. Setting up a basic 1-point mixture model and then running this on the object returns a model with the lowest AIC.

Progress will be printed to the screen.
}
\value{
  result: a \code{\link{ds.mixture}} model object of the best (AIC-wise) model. 
}
\references{ 
Miller, D.L. and L. Thomas (in prep.). Mixture model distance sampling detection functions.
 }
\author{David Lawrence Miller}
\seealso{\code{\link{summary.ds.mixture}}, \code{\link{step.ds.mixture}},\code{\link{plot.ds.mixture}}, \code{\link{sim.mix}}, \code{\link{mmds.gof}},\code{\link{ds.mixture}}, \code{\link{mmds.gof}}, \code{\link{fitmix}}, \code{\link{mmds.pars}}, \code{\link{step.ds.mixture}}}
\examples{

library(mmds)
set.seed(0)
## simulate some line transect data from a 2 point mixture
sim.dat<-sim.mix(c(-0.223,-1.897,inv.reparam.pi(0.3)),2,100,1)
## fit the model
fit.sim.dat.1<-fitmix(sim.dat,1,1)
## find best AIC model
step.ds.mixture(fit.sim.dat.1)
}
\keyword{Statistical Models, Distance Sampling}
