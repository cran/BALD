\name{degreesOfFreedom}
\alias{degreesOfFreedom}
\title{A generic function to plot and/or return the posterior of the degrees of freedom for the Student-t in  models in BALD.}
\description{A generic function to plot and/or return the posterior of the degrees of freedom for the Student-\eqn{t} in  models in \pkg{BALD}.}
\details{When there is zero skew, the degrees of freedom are the degrees of freedom for the non-skewed \eqn{t}.
See \code{vignette('BALD')}.}
\references{Kim, Y., and J. McCulloch (2007) \dQuote{The Skew-Student Distribution with Application to U.S. Stock Market Returns and the Equity Premium,} Department of Economics, Ohio State University, October 2007.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=degreesOfFreedom,AnnualAggLossDevModelOutput-method]{degreesOfFreedom("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the degrees of freedom.}
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
standard.model.output <- runLossDevModel(
standard.model.input,
burnIn=30.0E+3,
sampleSize=30.0E+3,
thin=10)
degreesOfFreedom(standard.model.output)
}}
