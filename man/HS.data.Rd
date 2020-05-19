\name{HS.data}
\alias{HS.data}
\docType{data}
\title{
Holizinger and Swineford's Data
}
\description{
This data set, for scores on a variety of tests, was originally in the MBESS package. A new version of the data set in that package doesn't appear to be identical to this one.
}
\usage{HS.data}
\format{
  A data frame with 301 observations on the following 32 variables.
  \describe{
    \item{\code{id}}{a numeric vector}
    \item{\code{Gender}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{grade}}{a numeric vector}
    \item{\code{agey}}{a numeric vector}
    \item{\code{agem}}{a numeric vector}
    \item{\code{school}}{a factor with levels \code{Grant-White} \code{Pasteur}}
    \item{\code{visual}}{a numeric vector}
    \item{\code{cubes}}{a numeric vector}
    \item{\code{paper}}{a numeric vector}
    \item{\code{flags}}{a numeric vector}
    \item{\code{general}}{a numeric vector}
    \item{\code{paragrap}}{a numeric vector}
    \item{\code{sentence}}{a numeric vector}
    \item{\code{wordc}}{a numeric vector}
    \item{\code{wordm}}{a numeric vector}
    \item{\code{addition}}{a numeric vector}
    \item{\code{code}}{a numeric vector}
    \item{\code{counting}}{a numeric vector}
    \item{\code{straight}}{a numeric vector}
    \item{\code{wordr}}{a numeric vector}
    \item{\code{numberr}}{a numeric vector}
    \item{\code{figurer}}{a numeric vector}
    \item{\code{object}}{a numeric vector}
    \item{\code{numberf}}{a numeric vector}
    \item{\code{figurew}}{a numeric vector}
    \item{\code{deduct}}{a numeric vector}
    \item{\code{numeric}}{a numeric vector}
    \item{\code{problemr}}{a numeric vector}
    \item{\code{series}}{a numeric vector}
    \item{\code{arithmet}}{a numeric vector}
    \item{\code{paperrev}}{a numeric vector}
    \item{\code{flagssub}}{a numeric vector}
  }
}

\source{
Originally from Holzinger and Swineford (1939). This copy is originally from version 4.6.0 of the MBESS package.
}

\references{
Holzinger, K. J. and Swineford, F. A. (1939). 
A study in factor analysis: The stability of a bi-factor solution. 
Supplementary Education Monographs, 48. University of Chicago.
}

\examples{
summary(HS.data)
}
\keyword{datasets}
