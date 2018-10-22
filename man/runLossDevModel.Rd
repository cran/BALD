\name{runLossDevModel}
\alias{runLossDevModel}
\title{A generic function to run models in BALD.}
\description{A generic function to run models in \pkg{BALD}.
See \code{vignette('BALD')}.
Overriding methods must return a valid output object.}
\value{object of class \code{LossDevModelOutput}.}
\seealso{\code{\link[=runLossDevModel,LossDevModelInput-method]{runLossDevModel("LossDevModelInput")}}}
\arguments{\item{object}{The object containing the model to estimate.}}
\examples{rm(list=ls())
library(BALD)
options(device.ask.default=FALSE)
library(BALD)
data(IncrementalGeneralLiablityTriangle)
IncrementalGeneralLiablityTriangle <- as.matrix(IncrementalGeneralLiablityTriangle)
print(IncrementalGeneralLiablityTriangle)
data(PCE)
PCE <- as.matrix(PCE)[,1]
PCE.rate <- PCE[-1] / PCE[-length(PCE)] - 1
PCE.rate.length <- length(PCE.rate)
print(PCE.rate[(-10):0 + PCE.rate.length])
PCE.years <- as.integer(names(PCE.rate))
years.available <- PCE.years <= max(as.integer(
dimnames(IncrementalGeneralLiablityTriangle)[[1]]))
PCE.rate <- PCE.rate[years.available]
PCE.rate.length <- length(PCE.rate)
print(PCE.rate[(-10):0 + PCE.rate.length])
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
}}
