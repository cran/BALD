\name{stochasticInflationRhoParameter}
\alias{stochasticInflationRhoParameter}
\title{A generic function to plot and/or return the posterior of the stochastic inflation rho parameter for models in BALD.}
\description{A generic function to plot and/or return the posterior of the stochastic inflation rho parameter for models in \pkg{BALD}.}
\details{If the model incorporates a stochastic rate of inflation, then that rate is assumed to follow (on the log scale) an autoregressive process of order 1.
(The autoregressive process of order 1 is the discrete equivalent to an Ornstein-Uhlenbeck process.)
This function plots the posterior for the \eqn{rho} parameter, assuming one was estimated.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=stochasticInflationRhoParameter,AnnualAggLossDevModelOutput-method]{stochasticInflationRhoParameter("AnnualAggLossDevModelOutput")}}
\code{\link{stochasticInflationStationaryMean}}
\code{\link{stochasticInflation}}}
\arguments{\item{object}{The object from which to plot and/or return the stochastic inflation \eqn{rho} parameter.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}}
\examples{rm(list=ls())
options(device.ask.default=FALSE)
library(BALD)
data(IncrementalGeneralLiablityTriangle)
IncrementalGeneralLiablityTriangle <- as.matrix(IncrementalGeneralLiablityTriangle)
print(IncrementalGeneralLiablityTriangle)
data(PCE)
PCE <- as.matrix(PCE)[,1]
PCE.rate <- PCE[-1] / PCE[-length(PCE)] - 1
PCE.rate.length <- length(PCE.rate)
PCE.years <- as.integer(names(PCE.rate))
years.available <- PCE.years <= max(as.integer(
dimnames(IncrementalGeneralLiablityTriangle)[[1]]))
PCE.rate <- PCE.rate[years.available]
PCE.rate.length <- length(PCE.rate)
standard.model.input <- makeStandardAnnualInput(
incremental.payments = IncrementalGeneralLiablityTriangle,
stoch.inflation.weight = 1,
non.stoch.inflation.weight = 0,
stoch.inflation.rate = PCE.rate,
exp.year.type = 'ay',
extra.dev.years=5,
use.skew.t=TRUE)
\dontrun{
standard.model.output <- runLossDevModel(
standard.model.input,
burnIn=30.0E+3,
sampleSize=30.0E+3,
thin=10)
stochasticInflationRhoParameter(standard.model.output)
}}
