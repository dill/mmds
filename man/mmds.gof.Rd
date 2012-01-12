\name{mmds.gof}
\alias{mmds.gof}
\title{Goodness of fit for mixture model detection functions}
\description{Goodness of fit testing for detection for mixture model detection functions. See Details below...
}
\details{
Two goodness of fit tests are provided: the Cramer-von Mises and the Kolmogorov-Smirnov. Both are implemented a in Buckland et al. (2004).

Print methods are provided, so accessing the \code{ks} and \code{cvm} elements of a \code{\link{ds.mixture}} object will give suitable summaries.
}
\references{ 
Miller, D.L. and L. Thomas (in prep.). Mixture model distance sampling detection functions.
 }
\author{David Lawrence Miller}
\seealso{\code{\link{summary.ds.mixture}}, \code{\link{step.ds.mixture}},\code{\link{plot.ds.mixture}}, \code{\link{sim.mix}}, \code{\link{mmds.gof}},\code{\link{ds.mixture}}, \code{\link{mmds.gof}}, \code{\link{fitmix}}, \code{\link{mmds.pars}}, \code{\link{step.ds.mixture}}}
\keyword{Statistical Models, Distance Sampling}
