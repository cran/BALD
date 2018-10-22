\name{triResi}
\alias{triResi}
\title{A generic function to plot and/or return residuals for models in the BALD package.}
\description{A generic function to plot and/or return residuals for models in the \pkg{BALD} package. See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=triResi,AnnualAggLossDevModelOutput-method]{triResi("AnnualAggLossDevModelOutput")}}
\code{\link{QQPlot}}}
\arguments{\item{object}{The object from which to plot and/or return the residuals.}
\item{standardize}{A logical value.  If \code{TRUE}, the plotted and returned residuals are normalized to their respective standard deviation.}
\item{timeAxis}{A character value describing along which of the three time axes to plot the residuals: \sQuote{dy} for development year time, \sQuote{cy} for calendar year time, \sQuote{ey} for exposure year time.}
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
#residule plot by Development Year
triResi(standard.model.output, timeAxis='dy')
#residule plot by exposure year
triResi(standard.model.output, timeAxis='ey')
#residule plot by calendar year
triResi(standard.model.output, timeAxis='cy')
}}
