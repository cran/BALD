\name{calendarYearEffectErrorTracePlot}
\alias{calendarYearEffectErrorTracePlot}
\title{A generic function to generate the trace plots for select calendar year effect errors.}
\description{A generic function to generate the trace plots for select calendar year effect errors.}
\details{The calendar year effect is comprised of two components: 1) a prior expected value that may be unique to every cell and 2) a diagonal-specific error term.
This function generates trace plots for the diagonal specific error terms only.
See \code{vignette('BALD')}.}
\value{NULL invisibly.  Only called for the side effect of plotting.}
\seealso{\code{\link[=calendarYearEffectErrorTracePlot,AnnualAggLossDevModelOutput-method]{calendarYearEffectErrorTracePlot("AnnualAggLossDevModelOutput")}}
\code{\link{calendarYearEffectErrors}}}
\arguments{\item{object}{The object from which to generate the trace plots.}
\item{elements}{A numeric vector indicating the elements for which to plot the trace.  Valid values are 2 through the total number of exposure years(observed and forecast).  If NULL, values are selected automatically.}}
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
calendarYearEffectErrorTracePlot(standard.model.output)
}}
