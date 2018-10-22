\name{accountForZeroPayments}
\alias{accountForZeroPayments}
\title{A function to take a triangle estimated without considering zero payments, and account for the possibility of zero payments.}
\usage{accountForZeroPayments(object, burnIn=1000, nAddapt=1000)}
\description{A function to take a triangle estimated without considering zero payments, and account for the possibility of zero payments.}
\details{As incremental payments are modeled on the log scale, zero payments (and negative payments) are treated as missing values.
So, without somehow accounting for zero payments, the estimated payments would be overstated.
Zero payments are accounted for by weighting the predicted payment (given that the payment is greater than zero) with the probability that this payment is zero.
(Negative payments are not (currently) accounted for.)
Currently the trajectory for this probably follows a gompertz curve and is constant across exposure years.
This is currently implemented as a function but may be switched to a method.
See \code{vignette('BALD')}.}
\arguments{\item{object}{The object containing the triangle estimated without accounting for zero payments.}
\item{burnIn}{An integer to represent the number of initial \acronym{MCMC} iterations to be discarded. (The adaptive phase (\code{nAddapt}) is not considered part of \code{burnIn}.)}
\item{nAddapt}{The length of the adaptive phase for the \acronym{MCMC} algorithm. (Default is \code{trunc(burnIn/4)+1}.)}}
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
}}
