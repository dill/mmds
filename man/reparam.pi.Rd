\name{reparam.pi}
\alias{reparam.pi}
\title{Reparametrize mixture proportions}
\usage{
reparam.pi(thetas)
}
\arguments{
  \item{thetas}{Mixture proportions in their
  parametrisation for optimization.}
}
\value{
a vector of parameters
}
\description{
Re-parameterise the mixture proportions so that when there
is more than a 2-point mixture, the proportions sum to 1.
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
inv.reparam.pi
}

