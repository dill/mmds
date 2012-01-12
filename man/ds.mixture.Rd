\name{ds.mixture}
\alias{ds.mixture}
\title{A fitted Mixture Model Detection Function Object}
\description{
The fitted mixture model detection function object returned by \code{\link{fitmix}}. Knowledge of most of this is not useful. Use \code{link{summary.ds.mixture}} for result summaries.}
\value{A \code{\link{ds.mixture}} object has the following elements:
   \item{distance}{Vector of distances used in the analysis.}
   \item{likelihood}{Value of the log-likelihood at the maxima.}
   \item{pars}{Parmeter estimates. See \code{\link{mmds.pars}} for more information.}
   \item{mix.terms}{Number of mixture terms fit.}
   \item{width}{Truncation distance used.}
   \item{z}{List containing the matrix of covariates used. Output from \code{\link{model.matrix}}.}
   \item{zdim}{Number of columns of \code{z}. See \code{\link{mmds.pars}} for more information.}
   \item{hessian}{Hessian matrix at the maxima.}
   \item{pt}{Logical indicating whether the data were from a point transect survey.}
   \item{data}{Data frame after truncation}
   \item{ftype}{Type of detection function.}
   \item{ctrl.options}{Options passed to \code{\link{optim}}.}
   \item{showit}{Debug level.}
   \item{opt.method}{Optimisation method used.}
   \item{usegrad}{Were analytic gradients used?}
   \item{model.formula}{Model formula.}
   \item{mu}{Per-observation effective trip width/effective area of detection.} 
   \item{pa.vec}{Vector of per-observation detectabilities.}
   \item{N}{Estimate of N in the covered area (Horvitz-Thompson).}
   \item{pa}{Average detectability.}
   \item{pars.se}{Standard errors of the parameters.}
   \item{N.se}{Standard error of the Horvitz-Thompson estimate of the abundance.}
   \item{pa.se}{Standard error of the average detectability.}
   \item{aic}{AIC of the fitted model.}
   \item{cvm}{Cramer-von Mises GoF test results. List containing: \code{p}, the p-value and \code{W}, the test statistic.}
   \item{ks}{Kolmogorov-Smirnov test results. List containing: \code{p}, the p-value and \code{Dn}, the test statistic. See \code{\link{mmds.gof}} for more information.}
}
\references{ 
Miller, D.L. and L. Thomas (in prep.). Mixture model distance sampling detection functions.
 }
\author{David Lawrence Miller}
\note{
\code{ds.mixture} objects can be passed to \code{\link{step.ds.mixture}} to select number of mixture components based on AIC score. 
}
\seealso{ \code{\link{summary.ds.mixture}}, \code{\link{step.ds.mixture}},\code{\link{plot.ds.mixture}}, \code{\link{sim.mix}}, \code{\link{mmds.gof}},\code{\link{ds.mixture}}, \code{\link{mmds.gof}}, \code{\link{mmds.pars}}, \code{\link{step.ds.mixture}}}
\keyword{Statistical Models, Distance Sampling}
