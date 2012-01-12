\name{plot.ds.mixture}
\alias{plot.ds.mixture}
\title{Plot the results of a mixture model detection function fit.}
\description{
Plots the detection function (or pdf) of a fitted mixture model detection function, optionally overlayed on a histogram of the observed data.
}
\usage{
\method{plot}{ds.mixture}(x,style="",main="",breaks="Sturges",ylim=NULL,xlim=NULL,pdf=FALSE,plot.formula=NULL,hide.hist=FALSE,nomf=FALSE,x.axis=NULL,...)
}
\arguments{
   \item{x}{a \code{\link{ds.mixture}} object.}
   \item{style}{If set to "comp", composite plots of the detection function will be shown for the detection function (the averaged detection function in the covariate case).}
   \item{main}{(A vector of) title(s) for the plot(s). By default these are set by the function (and are fairly ugly but descriptive).}
   \item{breaks}{Breaks to be used for the histogram. This can be a vector of numbers or any of the permissable options used in \code{\link{hist}}. Defaults to "Sturges".}
   \item{ylim}{Used to manually set the y limit of the plot. Defaults to NULL.}
   \item{xlim}{Used to manually set the x limit of the plot. Defaults to NULL.}
   \item{pdf}{Should the pdf be plotted rather than detection function be plotted? Only really useful with point transect data. Defaults to FALSE.}
   \item{plot.formula}{Formula of covariates to be plotted. Defaults to \code{NULL}, which plots all covariates. No effect with non-covariate models.}
   \item{hide.hist}{Should the histogram be hidden, leaving only the detection function (or pdf) to be plotted? Defaults to FALSE.}
   \item{nomf}{Should the mfrow value be altered? Useful when creating custom plots for publication. Defaults to FALSE (yes, change the mfrow value).}
   \item{x.axis}{Set the x axis labels. Again, useful for publication plots. Defaults to NULL, which uses the default R values.}
   \item{...}{Anything else!}
}
\details{
For covariate models, all the levels of factor variables are plotted or the 25, 50 and 75th percentiles of continuous variables are plotted averaged over the values of the other covariates.
}
\value{
  result: a plot. 
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
## plot 
plot.ds.mixture(fit.sim.dat) 
}
\keyword{Statistical Models, Distance Sampling}
