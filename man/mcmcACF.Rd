\name{mcmcACF}
\alias{mcmcACF}
\title{A generic function to plot autocorrelations found in the MCMC samples for select parameters.}
\description{A generic function to plot autocorrelations found in the \acronym{MCMC} samples for select parameters.}
\details{Chains with high autocorrelation require a longer burnin and more samples to fully explore the parameter space.
See \code{vignette('BALD')}.}
\value{Called for the side effect of plotting.}
\seealso{\code{\link[=mcmcACF,StandardAnnualAggLossDevModelOutput-method]{mcmcACF("StandardAnnualAggLossDevModelOutput")}}
\code{\link[=mcmcACF,BreakAnnualAggLossDevModelOutput-method]{mcmcACF("BreakAnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot autocorrelations.}}
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
mcmcACF(standard.model.output)
}}
