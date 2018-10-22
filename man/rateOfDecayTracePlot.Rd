\name{rateOfDecayTracePlot}
\alias{rateOfDecayTracePlot}
\title{A generic function to plot the trace plots for select rate of decay values.}
\description{A generic function to plot the trace plots for select rate of decay values.
See \code{vignette('BALD')}.}
\value{NULL invisibly.  Only called for the side effect of plotting.}
\seealso{\code{\link[=rateOfDecayTracePlot,StandardAnnualAggLossDevModelOutput-method]{rateOfDecayTracePlot("StandardAnnualAggLossDevModelOutput")}}
\code{\link{rateOfDecay}}}
\arguments{\item{object}{The object from which to generate the trace plots.}
\item{elements}{A numeric vector indicating for which elements to plot the trace.  Valid values are 2 through the number of columns in the observed triangle.  If NULL, values are selected automatically.}
\item{\dots}{Additional arguments used by methods.}}
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
rateOfDecayTracePlot(standard.model.output)
}}
