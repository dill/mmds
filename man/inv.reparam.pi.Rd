\name{inv.reparam.pi}
\alias{inv.reparam.pi}
\title{Reparametrize mixture proportions}
\usage{
  inv.reparam.pi(mix.prop, lastpar = FALSE)
}
\arguments{
  \item{mix.prop}{A set of mixture proportions.}

  \item{lastpar}{Is the last parameter provided, i.e. does
  sum(mix.prop)=1?}
}
\value{
  a vector of parameters
}
\description{
  Re-parameterise the mixture proportions so that when
  there is more than a 2-point mixture, the proportions sum
  to 1.
}
\section{Notes}{
  See Miller and Thomas for information on exactly how
  these are calculated. Thanks go to David Borchers for
  proposing the trick.
}
\examples{
library(mmds)
reparam.pi(inv.reparam.pi(0.3))
reparam.pi(inv.reparam.pi(c(0.3,0.4,0.1),lastpar=TRUE))
}
\author{
  David L. Miller
}
\references{
  Miller, D.L. and L. Thomas (in prep.). Mixture model
  distance sampling detection functions.
}
\seealso{
  reparam.pi
}

