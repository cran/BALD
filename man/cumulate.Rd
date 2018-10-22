\name{cumulate}
\alias{cumulate}
\title{A function to cumulate a triangle.}
\usage{cumulate(triangle)}
\description{A function to cumulate a triangle.}
\details{PC actuaries aggregate historical loss experience in a triangle two dimensions form.
The two dimensions are accident years and development years.
Most of the time, the loss experience is paid loss and/or incurred loss.
Incurred loss is the sum of paid loss and case reserve put up for the estimate of future payment of the claims.
For paid losses triangle, each cell could represent either net payment of the year or paid to date throughout the life of the claims.
The purpose of this function is to transform the incremental payment triangle to cumulative triangle.}
\value{A matrix resulting from cumulating the input triangle.
This is the cumulative triangles.
Each cell is the payment-to-date of the claims.}
\arguments{\item{triangle}{A matrix of incremental payments.
Or the incremental payment triangle.
Each cell of the incremental triangle is the payment during that specific year.}}
\examples{library(BALD) 
#load Generial Liability
data(IncrementalGeneralLiablityTriangle)
print(IncrementalGeneralLiablityTriangle)
#cumulate the triangle
cumulate(IncrementalGeneralLiablityTriangle)}
