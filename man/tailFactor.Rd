\name{tailFactor}
\alias{tailFactor}
\title{A generic function to plot and/or return the predicted tail factors for a specific attachment point.}
\description{A generic function to plot and/or return the predicted tail factors for a specific attachment point.}
\details{The tail factor is the ratio of the estimated ultimate loss to cumulative loss at some point in development time.
This is a generic function to allow for the retrieval and illustration the tail factor by exposure year.

\bold{Note on \code{firstIsHalfReport} and \code{attachment}:} \code{firstIsHalfReport} refers to the first column of the triangle.
For policy year triangles, the first column is often referred to as a \dQuote{half-report}, the second column is called \dQuote{first-report}, the third column is called \dQuote{second-report}, etc.
If \code{firstIsHalfReport=TRUE}, then \code{tailFactor} will assume the triangle is arranged in such a way that the first column is the \dQuote{half-report}
and \code{attachment=1} indicates that the charted tail factor attaches at the cumulative loss through the second column.  If \code{firstIsHalfReport=FALSE},
then \code{attachment=1} indicates that the charted tail factor attaches at the cumulative loss through the first column.  Since \code{attachment} must be coercible to an integer,
it is impossible to plot half-to-ultimate tail factors; however, they are the first column in the returned matrix.

\code{firstIsHalfReport} can be \code{NA} (the default)
if the exposure year type was specified to be one of \dQuote{policy year} or \dQuote{accident year} at the time the input object was constructed (see \code{\link{makeStandardAnnualInput}}
or \code{\link{makeBreakAnnualInput}}).  An exposure year type of \dQuote{policy year} corresponds to \code{firstIsHalfReport=TRUE},
and an exposure year type of \dQuote{accident year} corresponds to \code{firstIsHalfReport=FALSE}.  Setting \code{firstIsHalfReport} to a non-missing value will override this default.

If \code{expYearRange} is \dQuote{fullyObs}, then only exposure years with a non missing value in the first column will be plotted.
See \code{vignette('BALD')}.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=tailFactor,StandardAnnualAggLossDevModelOutput-method]{tailFactor("StandardAnnualAggLossDevModelOutput")}}
\code{\link[=tailFactor,BreakAnnualAggLossDevModelOutput-method]{tailFactor("BreakAnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot the predicted tail factors and return tail factors for \emph{all} attachment points.}
\item{attachment}{An integer value specifying the attachment point for the tail.  Must be at least 1. See Details for more information.}
\item{useObservedValues}{A logical value.  If \code{TRUE}, observed values are substituted for predicted values whenever possible in the calculation.  If \code{FALSE}, only predicted values are used.}
\item{firstIsHalfReport}{A logical value or \code{NA}.  See Details for more info.}
\item{finalAttachment}{An integer value must be at least 1. Default value is \code{attachment}.  A call to \code{tailFactor} will return (invisibly) a matrix of tail factors through this value.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}
\item{expYearRange}{Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.}}
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
tailFactor(standard.model.output,10)
}}
