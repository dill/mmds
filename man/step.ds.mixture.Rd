\name{step.ds.mixture}
\alias{step.ds.mixture}
\title{Stepwise selection of mixture components}
\usage{
  step.ds.mixture(ds.object, max.terms = 4)
}
\arguments{
  \item{ds.object}{\code{\link{ds.mixture}} object.}

  \item{max.terms}{Maximum number of mixture components to
  attempt to fit.  Default 4.}
}
\value{
  a \code{\link{ds.mixture}} model object of the best
  (AIC-wise) model.
}
\description{
  Uses AIC to select the number of mixture components.
}
\section{Details}{
  This routine is most useful during model building.
  Setting up a basic 1-point mixture model and then running
  this on the object returns a model with the lowest AIC.

  Progress will be printed to the screen.
}
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
\author{
  David L. Miller
}

