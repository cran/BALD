\name{standardDeviationVsDevelopmentTime}
\alias{standardDeviationVsDevelopmentTime}
\title{A generic function to plot and/or return the posterior estimated standard deviation by development year.}
\description{A generic function to plot and/or return the posterior estimated standard deviation by development year.}
\details{Aggregate loss development models in \pkg{BALD} allow for changes (by development year) in the measurement error around the log incremental payments.
This is a generic function that allows for the retrieval and illustration of this standard deviation.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=standardDeviationVsDevelopmentTime,AnnualAggLossDevModelOutput-method]{standardDeviationVsDevelopmentTime("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the estimated standard deviation by development year.}
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
standardDeviationVsDevelopmentTime(standard.model.output)
}}
