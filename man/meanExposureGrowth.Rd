\name{meanExposureGrowth}
\alias{meanExposureGrowth}
\title{A generic function to plot and/or return the posterior of the mean exposure growth for models in BALD.}
\description{A generic function to plot and/or return the posterior of the mean exposure growth for models in \pkg{BALD}.}
\details{(Optionally) exposure growth is modeled as an ar1 process.  This inherently assumes that periods of high exposure growth are (or at least have the possibility of being) followed by continued high periods.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=meanExposureGrowth,AnnualAggLossDevModelOutput-method]{meanExposureGrowth("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the mean exposure growth.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
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
standard.model.output <- runLossDevModel(standard.model.input,
burnIn=30.0E+3,
sampleSize=30.0E+3,
thin=10)
meanExposureGrowth(standard.model.output)
}}
