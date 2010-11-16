\name{ML.methods}
\alias{ML.methods}
\alias{GLS.methods}
\alias{anova.objectiveML}
\alias{logLik.objectiveML}
\alias{AIC.objectiveML}
\alias{AICc.objectiveML}
\alias{BIC.objectiveML}
\alias{CAIC.objectiveML}
\alias{deviance.objectiveML}
\alias{print.objectiveML}
\alias{print.objectiveGLS}
\alias{summary.objectiveML}
\alias{summary.objectiveGLS}
\alias{print.summary.objectiveML}

\title{
Methods for sem Objects Fit Using the \code{objectiveML} and \code{objectiveGLS} Objective Functions
}

\description{
These functions are for objects fit by \code{\link{sem}} using the \code{\link{objectiveML}} (multivariate-normal full-information maximum-likelihood) 
and \code{\link{objectiveGLS}} (generalized least squares) objective functions.
}

\usage{
\method{anova}{objectiveML}(object, model.2, robust=FALSE, ...)

\method{logLik}{objectiveML}(object, ...)
\method{deviance}{objectiveML}(object, ...)
\method{AIC}{objectiveML}(object, ..., k)
\method{AICc}{objectiveML}(object, ...)
\method{BIC}{objectiveML}(object, ...)
\method{CAIC}{objectiveML}(object, ...)

\method{print}{objectiveML}(x, ...)
\method{print}{objectiveGLS}(x, ...)
\method{summary}{objectiveML}(object, digits=5, conf.level=.90, robust=FALSE, analytic.se=object$t <= 100, ...)
\method{summary}{objectiveGLS}(object, ...)
}

\arguments{
  \item{object, model.2, x}{an object inheriting from class \code{objectiveML} or \code{objectiveGLS}.}
  \item{robust}{if \code{TRUE}, compute robust standard errors or test.}
  \item{k, \ldots}{ignored.}
  \item{digits}{digits to be printed.}
  \item{conf.level}{level for confidence interval for the RMSEA index (default is .9).}
  \item{analytic.se}{use analytic (as opposed to numeric) coefficient standard errors; default is \code{TRUE} is there are no more than
    100 parameters in the model and \code{FALSE} otherwise.}
}

\author{John Fox \email{jfox@mcmaster.ca} and Jarrett Byrnes}

\references{
See \code{\link{sem}}.
}

\seealso{
\code{\link{sem}}, \code{\link{objective.functions}}, \code{\link{modIndices.objectiveML}} 
}

\keyword{models}