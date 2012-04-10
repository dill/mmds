\name{sim.mix}
\alias{sim.mix}
\title{Simulate data from a mixture model detection function}
\usage{
  sim.mix(pars, mix.terms, n, width, zdim = 0, z = NULL,
    pt = FALSE, showit = FALSE)
}
\arguments{
  \item{pars}{Parameters of the model to fit. See
  \code{\link{mmds.pars}} for details.}

  \item{mix.terms}{Number of mixture components.}

  \item{n}{Number of data to generate.}

  \item{width}{Truncation distance.}

  \item{zdim}{Number of columns of \code{z}. Defaults to
  0.}

  \item{z}{Covariate data. Defaults to NULL. See details
  for more information.}

  \item{pt}{Should point transect data be generated?
  Defaults to FALSE.}

  \item{showit}{Print the acceptance rate. Defaults to
  FALSE.}
}
\value{
  a \code{data.frame} with the following columns:
  \tabular{ll}{ observed \tab Whether the object was
  observed, always \code{n} 1s. Kept for \code{mmds}
  compatability.\cr object \tab Object identifier, numbered
  1 to \code{n}. Kept for \code{mmds} compatability.\cr
  distance \tab Observed distances.\cr \tab Then follows as
  many columns as there are columns as \code{z}, named as
  in \code{z}.}
}
\description{
  Simulate data from a (line or point transect) mixture
  model detection function with or without covariates using
  rejection sampling.
}
\details{
  This routine uses rejection sampling, so may be rather
  slow of large sample sizes. Direct sampling will be
  available soon.
}
\examples{
library(mmds)
set.seed(0)
## simulate some line transect data from a 2 point mixture
sim.dat<-sim.mix(c(-0.223,-1.897,inv.reparam.pi(0.3)),2,100,1)
hist(sim.dat$distance)
}
\author{
  David L. Miller
}

