\name{gompertzParameters}
\alias{gompertzParameters}
\title{A generic function to plot and/or return the posterior of the parameters for the gompertz curve which describes the probability of payment.}
\description{A generic function to plot and/or return the posterior of the parameters for the gompertz curve which describes the probability of payment.}
\details{The scale parameter describes how steep the curve is.
Larger values are steeper.
Positive values indicate that the probability of a positive payment should decrease with development time.
(The scale is restricted to be positive.)

The fifty.fifty parameter gives the point (in development time) when the gompertz curve gives a probability of fifty percent.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link{gompertzParameters,AnnualAggLossDevModelOutputWithZeros-method}}}
\arguments{\item{object}{The object from which to plot and/or return the parameters.}
\item{parameter}{A character describing which parameter to plot. \dQuote{scale} for the scale parameter. \dQuote{fifty.fifty} for the point at which the gompertz give a probably of fifty percent.}
\item{plotDensity}{A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
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
break.model.output.w.zeros <- accountForZeroPayments(break.model.output)
gompertzParameters(break.model.output.w.zeros, parameter='scale')
}}
