\name{tsls}
\alias{tsls}
\alias{tsls.formula}
\alias{tsls.default}
\alias{fitted.tsls}
\alias{residuals.tsls}
\alias{coef.tsls}
\alias{vcov.tsls}
\alias{anova.tsls}
\alias{print.tsls}
\alias{summary.tsls}
\alias{print.summary.tsls}

\title{Two-Stage Least Squares}
\description{
  Fits a regression equation, such as an equation in a structural-equation model, by two-stage least squares. This
  is equivalent to direct instrumental-variables estimation when the number of
  instruments is equal to the number of predictors.
}
\usage{
\method{tsls}{formula}(formula, instruments, data, subset, weights, 
	na.action, contrasts=NULL, ...)
\method{tsls}{default}(y, X, Z, w, names=NULL, ...)

\method{print}{tsls}(x, ...)
\method{summary}{tsls}(object, digits=getOption("digits"), ...)
\method{print}{summary.tsls}(x, ...)
\method{anova}{tsls}(object, model.2, s2, dfe, ...)

\method{fitted}{tsls}(object, ...)
\method{residuals}{tsls}(object, ...)
\method{coef}{tsls}(object, ...)
\method{vcov}{tsls}(object, ...)
}

\arguments{
  \item{formula}{model formula for structural equation to be estimated; a regression
    constant is implied if not explicitly omitted.}
  \item{instruments}{one-sided model formula specifying instrumental variables.}
  \item{data}{an optional data frame containing the variables in the model. 
    By default the variables are taken from the environment from which \code{tsls} is
    called.}
  \item{subset}{an optional vector specifying a subset of observations to be used in 
    fitting the model.}
  \item{weights, w}{an optional vector of weights to be used in the fitting process; 
    if specified should be a non-negative numeric vector with one entry for each
    observation, to be used to compute weighted 2SLS estimates.}
  \item{na.action}{a function that indicates what should happen when the 
    data contain \code{NA}s. 
    The default is set by the \code{na.action} option.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} argument of 
    \code{\link[stats:model.matrix]{model.matrix.default}}.}
  \item{y}{Response-variable vector.}
  \item{X}{Matrix of predictors, including a constant (if one is in the model).}
  \item{Z}{Matrix of instrumental variables, including a constant (if one is in the model).}
  \item{names}{optional character vector of names for the columns of the \code{X} matrix.}
  \item{x, object, model.2}{objects of class \code{tsls} returned by \code{tsls.formula} 
    (or of class \code{summary.tsls}), for \code{anova} 
    containing nested models
    to be compared by an incremental \eqn{F}-test. One model should be nested in the other; the
    order of models is immaterial.}
  \item{s2}{an optional estimate of error variance for the denominator of the \eqn{F}-test. 
    If missing, the error-variance estimate is taken from the larger model.}
  \item{dfe}{optional error degrees of freedom, to be specified when an estimate of error
    variance is given.}
  \item{digits}{number of digits for summary output.}
  \item{...}{arguments to be passed down.}
}

\value{
  \code{tsls.formula} returns an object of class \code{tsls}, with the following components:
  \item{n}{number of observations.}
  \item{p}{number of parameters.}
  \item{coefficients}{parameter estimates.}
  \item{V}{estimated covariance matrix of coefficients.}
  \item{s}{residual standard error.}
  \item{residuals}{vector of residuals.}
  \item{response}{vector of response values.}
  \item{X}{model matrix.}
  \item{Z}{instrumental-variables matrix.}
  \item{response.name}{name of response variable, or expression evaluating to response.}
  \item{formula}{model formula.}
  \item{instruments}{one-sided formula for instrumental variables.}
}

\references{    
    Fox, J. (1979)
    Simultaneous equation models and two-stage least-squares.
    In Schuessler, K. F. (ed.)
    \emph{Sociological Methodology 1979}, Jossey-Bass.
    
    Greene, W. H. (1993)
    \emph{Econometric Analysis, Second Edition}, Macmillan.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{
summary(tsls(Q ~ P + D, ~ D + F + A, data=Kmenta))     # demand equation

summary(tsls(Q ~ P + F + A, ~ D + F + A, data=Kmenta)) # supply equation

anova(tsls(Q ~ P + F + A, ~ D + F + A, data=Kmenta),
    tsls(Q ~ 1, ~ D + F + A, data=Kmenta))
}
\keyword{models}
