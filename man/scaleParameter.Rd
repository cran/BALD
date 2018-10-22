\name{scaleParameter}
\alias{scaleParameter}
\title{A generic function to plot and/or return the posterior of the scale parameter for the Student-t measurement equation for models in BALD.}
\description{A generic function to plot and/or return the posterior of the scale parameter for the Student-\eqn{t} measurement equation for models in \pkg{BALD}.}
\details{As the degrees of freedom of the \eqn{t} goes to infinity, the scale parameter is the standard deviation of the resulting normal distribution (assuming zero skew).
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=scaleParameter,AnnualAggLossDevModelOutput-method]{scaleParameter("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the scale parameter.}
\item{column}{The scale parameter is allowed to vary with development time. Setting \code{column} results in the plotting and returning of the scale parameter corresponding to that column. Default value is \code{1}.}
\item{plotDensity}{A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
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
scaleParameter(standard.model.output)
}}
