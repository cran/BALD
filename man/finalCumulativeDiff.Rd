\name{finalCumulativeDiff}
\alias{finalCumulativeDiff}
\title{A generic function to plot and/or return the difference between final actual and predicted cumulative payments.}
\description{A generic function to plot and/or return the difference between final actual and predicted cumulative payments. See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=finalCumulativeDiff,AnnualAggLossDevModelOutput-method]{finalCumulativeDiff("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the difference.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}
\item{expYearRange}{Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.}}
\examples{rm(list=ls())
library(BALD)
options(device.ask.default=FALSE)
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
finalCumulativeDiff(standard.model.output)
}}
