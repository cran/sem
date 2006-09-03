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

\title{Two-Stage Least Squares}
\description{
  Fits an equation in a structural-equation model by two-stage least squares. This
  is equivalent to direct instrumental-variables estimation when the number of
  instruments is equal to the number of predictors.
}
\usage{
tsls(y, ...)
\method{tsls}{formula}(formula, instruments, data, subset, na.action, contrasts=NULL, ...)
\method{tsls}{default}(y, X, Z, names=NULL, ...)

\method{print}{tsls}(x, ...)
\method{summary}{tsls}(object, digits=4, ...)
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
    By default the variables are taken from the environment from which tsls is
    called.}
  \item{subset}{an optional vector specifying a subset of observations to be used in 
    fitting the model.}
  \item{na.action}{a function that indicates what should happen when the 
    data contain \code{NA}s. 
    The default is set by the \code{na.action} option.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of 
    \code{\link[stats:model.matrix]{model.matrix.default}}.}
  \item{y}{Response-variable vector.}
  \item{X}{Matrix of predictors, including a constant (if one is in the model).}
  \item{Z}{Matrix of instrumental variables, including a constant (if one is in the model).}
  \item{names}{optional character vector of names for the columns of the \code{X} matrix.}
  \item{x, object, model.2}{objects of class \code{tsls} returned by \code{tsls.formula}, 
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
data(Kmenta)
summary(tsls(Q ~ P + D, ~ D + F + A, data=Kmenta))     # demand equation

##  2SLS Estimates
##  
##  Model Formula: Q ~ P + D
##  
##  Instruments: ~D + F + A
##  
##  Residuals:
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -3.43e+00 -1.24e+00 -1.89e-01 -2.49e-13  1.58e+00  2.49e+00 
##  
##              Estimate Std. Error t value  Pr(>|t|)
##  (Intercept)  94.6333    7.92084  11.947 1.076e-09
##  P            -0.2436    0.09648  -2.524 2.183e-02
##  D             0.3140    0.04694   6.689 3.811e-06
##  
##  Residual standard error: 1.9663 on 17 degrees of freedom

summary(tsls(Q ~ P + F + A, ~ D + F + A, data=Kmenta)) # supply equation

##  2SLS Estimates
##  
##  Model Formula: Q ~ P + F + A
##  
##  Instruments: ~D + F + A
##  
##  Residuals:
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -4.87e+00 -1.26e+00  6.42e-01 -5.64e-12  1.47e+00  3.49e+00 
##  
##              Estimate Std. Error t value  Pr(>|t|)
##  (Intercept)  49.5324   12.01053   4.124 7.954e-04
##  P             0.2401    0.09993   2.402 2.878e-02
##  F             0.2556    0.04725   5.410 5.785e-05
##  A             0.2529    0.09966   2.538 2.193e-02
##  
##  Residual standard error: 2.4576 on 16 degrees of freedom


anova(tsls(Q ~ P + F + A, ~ D + F + A, data=Kmenta),
    tsls(Q ~ 1, ~ D + F + A, data=Kmenta))

##  Analysis of Variance
##  
##  Model 1:  Q ~ P + F + A   Instruments: ~D + F + A
##  Model 2:  Q ~ 1   Instruments: ~D + F + A
##  
##          Res.Df     RSS Df Sum of Sq      F  Pr(>F)
##  Model 1     16  96.633                            
##  Model 2     19 268.114  3   171.481 4.0507 0.02553

}
\keyword{models}
