\name{Tests}
\alias{Tests}
\docType{data}
\title{
Six Mental Tests
}
\description{
These data are from the SAS manual and consist of six mental tests for 32 students,
with some missing data. The three \code{x} variables are intended to load on a verbal
factor, and the three \code{y} variables on a math factor. The data can be used to
illustrate the estimation of a confirmatory factor analysis model by multinormal full-information
maximum-likelihood in the presence of missing data.
}

\usage{Tests}
\format{
  A data frame with 32 observations on the following 6 variables.
  \describe{
    \item{\code{x1}}{score on verbal test 1.}
    \item{\code{x2}}{score on verbal test 2.}
    \item{\code{x3}}{score on verbal test 3.}
    \item{\code{y1}}{score on math test 1.}
    \item{\code{y2}}{score on math test 2.}
    \item{\code{y3}}{score on math test 3.}
  }
}

\source{
Example 25.13 from \emph{SAS/STAT 9.22 User's Guide}, SAS Institute, 2010.
}

\keyword{datasets}
