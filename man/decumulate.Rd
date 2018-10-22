\name{decumulate}
\alias{decumulate}
\title{A function to decumulate a triangle.}
\usage{decumulate(triangle)}
\description{A function to decumulate a triangle.}
\value{A matrix resulting from decumulating the input triangle.}
\arguments{\item{triangle}{A matrix of cumulative payments.}}
\examples{library(BALD) 
#load Generial Liability
data(IncrementalGeneralLiablityTriangle)
#cumulate the incremental triangle
CumulateGeneralLiablityTriangle<-cumulate(IncrementalGeneralLiablityTriangle)
#decumulate the loss triangle
decumulate(CumulateGeneralLiablityTriangle)}
