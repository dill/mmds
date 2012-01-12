\name{fitmix}
\alias{fitmix}
\title{ Mixture Model Distance Detection Function Fitting}
\description{
Fits a mixture of half-normals as a detection function to distance sampling data collected via either line or point transects, possibly with covariates.
}
\usage{
fitmix(data,width,mix.terms=1,pt=FALSE,model.formula="~1",initialvalues=NULL,showit=0,ctrl.options=c(maxit=10000),opt.method="BFGS+SANN",usegrad=TRUE,ftype="hn")
}
\arguments{
   \item{data}{\code{data.frame} containing the distances and covariates to be used in the analusis (see Details).}
   \item{width}{Truncation distance.}
   \item{mix.terms}{Number of mixture components to use. Defaults to 1 (ie. CDS).}
   \item{pt}{Is the data from point transects? Default FALSE.}
   \item{model.formula}{Formula to be used for the covariates. Defaults to "~1" (ie. no covariates).}
   \item{initialvalues}{User supplied initialvalues if needed. Defaults to \code{NULL}. See \code{\link{mmds.pars}} for more information.}
   \item{showit}{Debugging level from 0 to 3, with 3 being most verbose. Defaults to 0.}
   \item{ctrl.options}{Options to give to the \code{\link{optim}}. Defaults to c(maxit=10000).}
   \item{opt.method}{Optimisation method to use, one of "BFGS", "BFGS+SANN" or "EM". Defaults to "BFGS+SANN", see Details.}
   \item{usegrad}{Should analytic derivatives be used in the optimisation? Default TRUE.}
   \item{ftype}{Function type to be used as the detection function, currently only "hn" only.}
}
\details{
This is the main routine that fits mixture model detection functions.

\code{data} should be a \code{data.frame} with (at least) a column named \code{distance}. Any covariates given in \code{model.formula} should be named in \code{data}. Note that rows with \code{distance} greater than \code{width} will be discarded.

See \code{\link{step.ds.mixture}} for AIC selection for the number of mixture components.
}
\value{
  result: a \code{\link{ds.mixture}} model object. 
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
fit.sim.dat<-fitmix(sim.dat,1,2)
## what happened?
summary(fit.sim.dat) 
}
\keyword{Statistical Models, Distance Sampling}
