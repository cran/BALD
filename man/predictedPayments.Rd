\name{predictedPayments}
\alias{predictedPayments}
\title{A generic function to plot predicted vs actual payments for models from the BALD package.}
\description{A generic function to plot predicted vs actual payments for models from the \pkg{BALD} package.}
\details{Because the model is Bayesian, each estimated payment comes as a distribution.
The median of this distribution is used as a point estimate when plotting and/or returning values.
Note: One cannot calculate the estimated incremental payments from the estimated cumulative payments (and vice versa) since the median of sums need not be equal to the sum of medians.

If \code{mergePredictedWithObserved=TRUE} and \code{type="incremental"}, then any observed incremental payment will be used in place of its corresponding incremental payment.
If \code{mergePredictedWithObserved=TRUE} and \code{type="cumulative"}, then only predicted incremental payments (by row) to the right of the last observed cumulative value will enter the calculation.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=predictedPayments,AnnualAggLossDevModelOutput-method]{predictedPayments("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot predicted vs actual payments and from which to return predicted payments.}
\item{type}{A single character value specifying whether to plot/return the predicted incremental or cumulative payments. Valid values are \dQuote{incremental} or \dQuote{cumulative.}  See details as to why these may not match up.}
\item{logScale}{A logical value.  If \code{TRUE}, then values are plotted on a log scale.}
\item{mergePredictedWithObserved}{A logical value.  See details.}
\item{plotObservedValues}{A logical value.  If \code{FALSE}, then only the predicted values are plotted.}
\item{plotPredictedOnlyWhereObserved}{A logical value.  If \code{TRUE}, then only the predicted incremental payments with valid corresponding observed (log) incremental payment are plotted. Ignored for \code{type="cumulative"}.}
\item{quantiles}{A vector of quantiles for the predicted payments to return.  Useful for constructing credible intervals.}
\item{plot}{A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
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
predictedPayments(standard.model.output)
}}
