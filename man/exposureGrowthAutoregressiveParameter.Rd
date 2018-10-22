\name{exposureGrowthAutoregressiveParameter}
\alias{exposureGrowthAutoregressiveParameter}
\title{A generic function to plot and/or return the posterior of the autoregressive parameter for the exposure growth for models in BALD.}
\description{A generic function to plot and/or return the posterior of the autoregressive parameter for the exposure growth for models in \pkg{BALD}.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=exposureGrowthAutoregressiveParameter,AnnualAggLossDevModelOutput-method]{exposureGrowthAutoregressiveParameter("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the autoregressive parameter which is associated with exposure growth.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
\examples{rm(list=ls())
library(BALD)
data(CumulativeAutoBodilyInjuryTriangle)
CumulativeAutoBodilyInjuryTriangle <- as.matrix(CumulativeAutoBodilyInjuryTriangle)
sample.col <- (dim(CumulativeAutoBodilyInjuryTriangle)[2] - 6:0)
print(decumulate(CumulativeAutoBodilyInjuryTriangle)[1:7, sample.col])
data(HPCE)
HPCE <- as.matrix(HPCE)[,1]
HPCE.rate <- HPCE[-1] / HPCE[-length(HPCE)] - 1
print(HPCE.rate[(-10):0 + length(HPCE.rate)])
HPCE.years <- as.integer(names(HPCE.rate))
max.exp.year <- max(as.integer(dimnames(CumulativeAutoBodilyInjuryTriangle)[[1]]))
years.to.keep <- HPCE.years <=  max.exp.year + 3
HPCE.rate <- HPCE.rate[years.to.keep]
break.model.input.w.ar1 <- makeBreakAnnualInput(
cumulative.payments = CumulativeAutoBodilyInjuryTriangle,
stoch.inflation.weight = 1,
non.stoch.inflation.weight = 0,
stoch.inflation.rate = HPCE.rate,
first.year.in.new.regime = c(1986, 1987),
prior.for.first.year.in.new.regime=c(2,1),
exp.year.type = 'ay',
extra.dev.years = 5,
use.skew.t = TRUE,
bound.for.skewness.parameter=5,
use.ar1.in.exposure.growth = TRUE)
\dontrun{
break.model.output.w.ar1 <- runLossDevModel(
break.model.input.w.ar1,
burnIn=30.0E+3,
sampleSize=30.0E+3,
thin=10)
exposureGrowthAutoregressiveParameter(break.model.output.w.ar1)
}}
