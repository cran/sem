\name{raw.moments}
\alias{raw.moments}
\alias{raw.moments.formula}
\alias{raw.moments.default}
\alias{print.rawmoments}

\title{Compute Raw Moments Matrix}
\description{
 Computes the ``uncorrected'' sum-of-squares-and-products matrix divided by the
 number of observations.
}
\usage{
raw.moments(object, ...)

\method{raw.moments}{formula}(formula, data, subset, na.action, 
    contrasts = NULL, ...)

\method{raw.moments}{default}(object, ...)

\method{print}{rawmoments}(x, ...)
}

\arguments{
  \item{object}{a one-sided model formula or an object coercible to a
    numeric matrix.}
  \item{formula}{a one-sided model formula specifying the model matrix for
    which raw moments are to be computed; note that a constant is included
    if it is not explicitly suppressed by putting \code{-1} in the formula.}
  \item{data}{an optional data frame containing the variables in the formula. 
    By default the variables are taken from the environment from which 
    \code{raw.moments} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be 
    used in computing moments.}
  \item{na.action}{a function that indicates what should happen when the data 
    contain \code{NA}s. The default is set by the \code{na.action} option.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of 
    \code{\link[stats:model.matrix]{model.matrix.default}}}.
  \item{x}{an object of class \code{rawmoments} to print.}
  \item{\dots}{arguments passed down.}
}

\value{
  Returns an object of class \code{rawmoments}, which is simply a matrix
  with an attribute \code{"N"} that contains the number of observations on
  which the moments are based.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{
data(Kmenta)
raw.moments(cbind(1, Kmenta))
raw.moments(~ Q + P + D + F + A, data=Kmenta)
}

\keyword{manip}

