\name{summary.ds.mixture}
\alias{summary.ds.mixture}
\title{Summarize a ds.mixture object}
\usage{
  \method{summary}{ds.mixture} (object, ...)
}
\arguments{
  \item{object}{A fitted mixture model detection function
  object.}

  \item{...}{Anything, but it will be ignored.}
}
\value{
  a summary of a \code{\link{ds.mixture}} object.
}
\description{
  Summarize a ds.mixture object. The function provides
  information on parameter estimates, estimates of the
  abundance in the covered area and the average
  detectability and their respective standard errors and
  coefficients of variation.
}
\examples{
library(mmds)
set.seed(0)
## simulate some line transect data from a 2 point mixture
sim.dat<-sim.mix(c(-0.223,-1.897,inv.reparam.pi(0.3)),2,100,1)
## fit the model
fit.sim.dat<-fitmix(sim.dat,1,2)
## what happened?
summary(fit.sim.dat)
}
\author{
  David L. Miller
}
\references{
  Miller, D.L. and L. Thomas (in prep.). Mixture model
  distance sampling detection functions.
}

