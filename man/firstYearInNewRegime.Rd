\name{firstYearInNewRegime}
\alias{firstYearInNewRegime}
\title{A generic function to plot and/or return the posterior change point.}
\description{A generic function to plot and/or return the posterior change point.}
\details{When incorporating a structural break, the user has the option of specifying either 1) the first year in which the new regime applies or 2) a (inclusive) range in which the first year in the new regime applies.
If the user specifies a range, the actual year is estimated as a model parameter.
This function allows for the retrieval/illustration of the posterior for this estimate.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=firstYearInNewRegime,BreakAnnualAggLossDevModelOutput-method]{firstYearInNewRegime("BreakAnnualAggLossDevModelOutput")}}
\code{\link{firstYearInNewRegimeTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the posterior change point estimate.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
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
firstYearInNewRegime(break.model.output)
}}
