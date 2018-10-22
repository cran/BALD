\name{firstYearInNewRegimeTracePlot}
\alias{firstYearInNewRegimeTracePlot}
\title{A generic function to generate the trace plot for the posterior change point.}
\description{A generic function to generate the trace plot for the posterior change point.
See \code{vignette('BALD')}.}
\value{Only called for the side effect of plotting.}
\seealso{\code{\link[=firstYearInNewRegimeTracePlot,BreakAnnualAggLossDevModelOutput-method]{firstYearInNewRegimeTracePlot("BreakAnnualAggLossDevModelOutput")}}
\code{\link{firstYearInNewRegimeTracePlot}}}
\arguments{\item{object}{The object from which to generate the trace plot for the change point estimate.}}
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
max.exp.year <- max(as.integer(
dimnames(CumulativeAutoBodilyInjuryTriangle)[[1]]))
years.to.keep <- HPCE.years <=  max.exp.year + 3
HPCE.rate <- HPCE.rate[years.to.keep]
break.model.input <- makeBreakAnnualInput(
cumulative.payments = CumulativeAutoBodilyInjuryTriangle,
stoch.inflation.weight = 1,
non.stoch.inflation.weight = 0,
stoch.inflation.rate = HPCE.rate,
first.year.in.new.regime = c(1986, 1987),
prior.for.first.year.in.new.regime=c(2,1),
exp.year.type = 'ay',
extra.dev.years = 5,
use.skew.t = TRUE,
bound.for.skewness.parameter=5)
\dontrun{
break.model.output <- runLossDevModel(
break.model.input,
burnIn=30.0E+3,
sampleSize=30.0E+3,
thin=10)
firstYearInNewRegimeTracePlot(break.model.output)
}}
