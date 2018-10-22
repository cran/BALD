\name{numberOfKnots}
\alias{numberOfKnots}
\title{A generic function to plot and/or return the posterior number of knots.}
\description{A generic function to plot and/or return the posterior number of knots.}
\details{The \link[=consumptionPath]{consumption path} (or calendar year effect and exposure growth adjusted log incremental payments) is modeled as a linear spline.
The number of knots (or places where the spline changes slope) in this spline is endogenous to the model and estimated by way of Reversible Jump Markov Chain Monte Carlo simulation.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.  Also returns statics on the number of knots.  Returned invisibly.}
\seealso{\code{\link{consumptionPath}}
\code{\link[=numberOfKnots,StandardAnnualAggLossDevModelOutput-method]{numberOfKnots("StandardAnnualAggLossDevModelOutput")}}
\code{\link[=numberOfKnots,BreakAnnualAggLossDevModelOutput-method]{numberOfKnots("BreakAnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot the number of knots.}
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
numberOfKnots(standard.model.output,10)
}}
