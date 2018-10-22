\name{calendarYearEffect}
\alias{calendarYearEffect}
\title{A generic function to plot and/or return the predicted and forecast calendar year effects for models in BALD.}
\description{A generic function to plot and/or return the predicted and forecast calendar year effects for models in \pkg{BALD}.}
\details{The calendar year effect is comprised of two components: 1) a prior expected value that may be unique to every cell (subject to weights and bounds) and 2) a diagonal-specific error term.
This function plots and returns the factor resulting from the combined effect of these two, which includes an autoregressive component if the model is estimated with such a feature.

The first cell is \code{NA}. Values in the first column represent the rate of inflation/escalation to the corresponding cell from the cell in the same column but previous row.
Values in the 2nd column and beyond represent the rate of inflation/escalation to the corresponding cell from the cell in the same row but previous column.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=calendarYearEffect,AnnualAggLossDevModelOutput-method]{calendarYearEffect("AnnualAggLossDevModelOutput")}}
\code{\link{calendarYearEffectErrors}}
\code{\link{autoregressiveParameter}}
\code{\link{standardDeviationOfCalendarYearEffect}}
\code{\link{calendarYearEffectErrorTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the calendar year effect.}
\item{restrictedSize}{A logical value.  If \code{TRUE}, the plotted calendar year effect is restricted to the square of dimension equal to the observed triangle with which the model was estimated.}
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
calendarYearEffect(standard.model.output)
}}
