\name{lossDevelopmentFactors}
\alias{lossDevelopmentFactors}
\title{A generic function to plot and/or return a table of predicted age-to-age loss development factors (or link ratios).}
\description{A generic function to plot and/or return a table of predicted age-to-age loss development factors (or link ratios).}
\details{While the model estimates ultimate losses directly, comparisons of predicted to observed development factors can give the user a better feel for the model's adequacy.
Since the model is Bayesian, each development factor comes as a distribution.  Only the median, as a point estimate, are plotted/returned.

The age-to-age factors are the ratios of the cumulative paid values at one period to the previous period.
Note that the median of products is not the product of medians, and thus it is possible (or rather likely) that age-to-age factors will not line up with age-to-ultimate factors (see \code{\link{tailFactor}}).
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.  Also returns a numeric matrix of plotted statistics.}
\seealso{\code{\link[=lossDevelopmentFactors,AnnualAggLossDevModelOutput-method]{lossDevelopmentFactors("AnnualAggLossDevModelOutput")}}
\code{\link{tailFactor}}}
\arguments{\item{object}{The object from which to plot and/or return loss development factors.}
\item{cex.text}{The \code{cex} value supplied to \code{text}. Adjusts the relative size of text.}
\item{linespace}{Adjusts the spacing between observed and predicted values.}
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
lossDevelopmentFactors(standard.model.output)
}}
