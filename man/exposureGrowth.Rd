\name{exposureGrowth}
\alias{exposureGrowth}
\title{A generic function to plot and/or return the posterior predicted exposure growth (corresponding to eta in the model).}
\description{A generic function to plot and/or return the posterior predicted exposure growth (corresponding to \emph{eta} in the model). See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting the exposure growth.  Also returns a named numeric vector for the median of the posterior for the exposure growth on the real (not log) scale.  Returned invisibly.}
\seealso{\code{\link[=exposureGrowth,AnnualAggLossDevModelOutput-method]{exposureGrowth("AnnualAggLossDevModelOutput")}}
\code{\link{exposureGrowthTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the posterior predicted exposure growth.}
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
years.available <- PCE.years <= max(
as.integer(dimnames(IncrementalGeneralLiablityTriangle)[[1]]))
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
exposureGrowthTracePlot(standard.model.output)
}}
