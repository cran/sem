\name{rawMoments}
\alias{rawMoments}
\alias{rawMoments.formula}
\alias{rawMoments.default}
\alias{cov2raw}
\alias{print.rawmoments}

\title{Compute Raw Moments Matrix}
\description{
 Computes the \dQuote{uncorrected} sum-of-squares-and-products matrix divided by the
 number of observations.
}
\usage{
\method{rawMoments}{formula}(formula, data, subset, na.action, 
    contrasts=NULL, ...)

\method{rawMoments}{default}(object, na.rm=FALSE, ...)

cov2raw(cov, mean, N, sd)

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
    \code{rawMoments} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be 
    used in computing moments.}
  \item{na.action}{a function that indicates what should happen when the data 
    contain \code{NA}s. The default is set by the \code{na.action} option.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} argument of 
    \code{\link[stats:model.matrix]{model.matrix.default}}}.
  \item{na.rm}{if \code{TRUE}, any data rows with missing data will be removed.}
  \item{cov}{a covariance or correlation matrix.}
  \item{mean}{a vector of means.}
  \item{N}{the number of observations on which the covariances or correlations 
    are based.}
  \item{sd}{an optional vector of standard deviations, to be given if \code{cov} is a
    correlation matrix.}
  \item{x}{an object of class \code{rawmoments} to print.}
  \item{\dots}{arguments passed down.}
}

\value{
  \code{rawMoments} and \code{cov2raw} return an object of class \code{rawmoments}, 
  which is simply a matrix
  with an attribute \code{"N"} that contains the number of observations on
  which the moments are based.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{
# the following are all equivalent (with the exception of the name of the intercept):

rawMoments(cbind(1, Kmenta))

rawMoments(~ Q + P + D + F + A, data=Kmenta)

Cov <- with(Kmenta, cov(cbind(Q, P, D, F, A)))
cov2raw(Cov, colMeans(Kmenta), nrow(Kmenta))
}

\keyword{manip}

