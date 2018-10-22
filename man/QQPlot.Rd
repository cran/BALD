\name{QQPlot}
\alias{QQPlot}
\title{A generic function to plot a Q-Q plot for models in the BALD package.}
\description{A generic function to plot a Q-Q plot for models in the \pkg{BALD} package.}
\details{This function plots sorted observed log incremental payments vs sorted predicted log incremental payments.
Credible intervals are also plotted. See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=QQPlot,AnnualAggLossDevModelOutput-method]{QQPlot("AnnualAggLossDevModelOutput")}}
\code{\link{triResi}}}
\arguments{\item{object}{The object from which to plot the values.}}
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
QQPlot(standard.model.output)
}}
