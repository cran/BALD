\name{stochasticInflation}
\alias{stochasticInflation}
\title{A generic function to plot and/or return predicted and forecast stochastic inflation rates for models in BALD.}
\description{A generic function to plot and/or return predicted and forecast stochastic inflation rates for models in \pkg{BALD}.}
\details{If the model incorporates a stochastic rate of inflation, then that rate is assumed to follow (on the log scale) an autoregressive process of order 1.
(The autoregressive process of order 1 is the discrete equivalent to an Ornstein-Uhlenbeck process.)
This function plots the median of the posterior predictive distribution for stochastic inflation (not on the log scale) rates by year.
Values are returned prior to the application of any limits or weights.
Note that for years where observed values are supplied, the model takes those values at face value.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=stochasticInflation,AnnualAggLossDevModelOutput-method]{stochasticInflation("AnnualAggLossDevModelOutput")}}
\code{\link{stochasticInflationRhoParameter}}
\code{\link{stochasticInflationStationaryMean}}}
\arguments{\item{object}{The object from which to plot and/or return the stochastic inflation rates.}
\item{extraYears}{An integer expressing the (maximum) number of years to plot (beyond the final observed year).  Must be at least zero.  Default is 15.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
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
stochasticInflation(standard.model.output)
}}
