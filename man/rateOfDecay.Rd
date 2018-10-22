\name{rateOfDecay}
\alias{rateOfDecay}
\title{A generic function to plot and/or return the esimtated rate of decay vs development year time.}
\description{A generic function to plot and/or return the esimtated rate of decay vs development year time.}
\details{The simplest definition of the rate of decay is the exponentiated first difference of the \link[=consumptionPath]{consumption path}.
This is a generic function to allow for the retrieval and illustration of the rate of decay.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.}
\seealso{\code{\link[=rateOfDecay,StandardAnnualAggLossDevModelOutput-method]{rateOfDecay("StandardAnnualAggLossDevModelOutput")}}
\code{\link[=rateOfDecay,BreakAnnualAggLossDevModelOutput-method]{rateOfDecay("BreakAnnualAggLossDevModelOutput")}}
\code{\link{consumptionPath}}
\code{\link{rateOfDecayTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the estimated rate of decay.}
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
rateOfDecay(standard.model.output)
}}
