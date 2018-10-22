\name{consumptionPath}
\alias{consumptionPath}
\title{A generic function to plot and/or return the estimated consumption path vs development year time.}
\description{A generic function to plot and/or return the estimated consumption path vs development year time.}
\details{At the heart of aggregate loss development models in \pkg{BALD} is the consumption path.
The consumption path is (on a log scale) the trajectory of incremental payments absent calendar year effects and with exposure normalized to the first row.
Note that the measurement error term is (possibly) a skewed \eqn{t} and as such (possibly) has a non zero mean.   The consumption path is absent any such shifts due to skewness.
This is a generic function that allows for the retrieval and illustration of this consumption path.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.}
\seealso{\code{\link[=consumptionPath,StandardAnnualAggLossDevModelOutput-method]{consumptionPath("StandardAnnualAggLossDevModelOutput")}}
\code{\link[=consumptionPath,BreakAnnualAggLossDevModelOutput-method]{consumptionPath("BreakAnnualAggLossDevModelOutput")}}
\code{\link{consumptionPathTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the estimated consumption path.}
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
consumptionPath(standard.model.output)
}}
