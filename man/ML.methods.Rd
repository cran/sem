\name{ML.methods}
\alias{ML.methods}
\alias{GLS.methods}
\alias{anova.objectiveML}
\alias{anova.objectiveFIML}
\alias{logLik.objectiveML}
\alias{logLik.objectiveFIML}
\alias{AIC.objectiveML}
\alias{AIC.objectiveFIML}
\alias{AICc.objectiveML}
\alias{AICc.objectiveFIML}
\alias{BIC.objectiveML}
\alias{BIC.objectiveFIML}
\alias{CAIC.objectiveML}
\alias{CAIC.objectiveFIML}
\alias{deviance.objectiveML}
\alias{deviance.objectiveFIML}
\alias{print.objectiveML}
\alias{print.objectiveGLS}
\alias{print.objectiveFIML}
\alias{summary.objectiveML}
\alias{summary.objectiveGLS}
\alias{summary.objectiveFIML}
\alias{print.summary.objectiveML}
\alias{print.msemObjectiveML}
\alias{print.msemObjectiveGLS}
\alias{summary.msemObjectiveML}
\alias{summary.msemObjectiveGLS}
\alias{deviance.msemObjectiveML}
\alias{AIC.msemObjectiveML}
\alias{AICc.msemObjectiveML}
\alias{BIC.msemObjectiveML}
\alias{deviance.msemObjectiveML}
\alias{anova.msemObjectiveML}
\alias{logLik.msemObjectiveML}


\title{
Methods for sem Objects Fit Using the \code{objectiveML}, \code{objectiveGLS}, \code{objectiveFIML}, \code{msemObjectiveML},
and \code{msemObjectiveGLS} Objective Functions
}

\description{
These functions are for objects fit by \code{\link{sem}} using the \code{\link{objectiveML}} (multivariate-normal full-information maximum-likelihood), \code{link{objectiveFIML}} (multivariate-normal full-information maximum-likihood in
the presence of missing data),
\code{\link{objectiveGLS}} (generalized least squares), and \code{\link{msemObjectiveML}} (multigroup multivariate-normal FIML) objective functions.
}

\usage{
\method{anova}{objectiveML}(object, model.2, robust=FALSE, ...)
\method{anova}{objectiveFIML}(object, model.2, ...)

\method{logLik}{objectiveML}(object, ...)
\method{logLik}{objectiveFIML}(object, saturated=FALSE, 
    intercept="Intercept", iterlim=1000, ...)
\method{deviance}{objectiveML}(object, ...)
\method{deviance}{objectiveFIML}(object, saturated.logLik, ...)
\method{deviance}{msemObjectiveML}(object, ...) 
\method{AIC}{objectiveML}(object, ..., k)
\method{AIC}{objectiveFIML}(object, saturated.logLik, ..., k)
\method{AIC}{msemObjectiveML}(object, ..., k)
\method{AICc}{objectiveML}(object, ...)
\method{AICc}{objectiveFIML}(object, saturated.logLik, ...)
\method{AICc}{msemObjectiveML}(object, ...)
\method{BIC}{objectiveML}(object, ...)
\method{BIC}{objectiveFIML}(object, saturated.logLik, ...)
\method{BIC}{msemObjectiveML}(object, ...)
\method{CAIC}{objectiveML}(object, ...)
\method{CAIC}{objectiveFIML}(object, saturated.logLik, ...)

\method{print}{objectiveML}(x, ...)
\method{print}{objectiveGLS}(x, ...)
\method{print}{objectiveFIML}(x, saturated=FALSE, ...)
\method{print}{msemObjectiveML}(x, ...)
\method{print}{msemObjectiveGLS}(x, ...)

\method{summary}{objectiveML}(object, digits=getOption("digits"), 
    conf.level=.90, robust=FALSE, analytic.se=object$t <= 500, 
    fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI",
      "IFI", "SRMR", "AIC", "AICc", "BIC", "CAIC"), ...)
\method{summary}{objectiveFIML}(object, digits=getOption("digits"), conf.level=.90,
    fit.indices=c("AIC", "AICc", "BIC", "CAIC"),
    saturated=FALSE, intercept="Intercept", saturated.logLik, ...)
\method{summary}{objectiveGLS}(object, digits=getOption("digits"), conf.level=.90, 
    fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR"),
    ...)
\method{summary}{msemObjectiveML}(object, digits=getOption("digits"), 
    conf.level=.90, robust=FALSE, 
    analytic.se=object$t <= 500,
    fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", 
      "IFI", "SRMR", "AIC", "AICc", "BIC"), ...)
\method{summary}{msemObjectiveGLS}(object, digits=getOption("digits"), 
    conf.level=.90,
    fit.indices=c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", 
      "CFI", "RNI", "IFI", "SRMR"), ...)
}

\arguments{
  \item{object, model.2, x}{an object inheriting from class \code{objectiveML}, \code{objectiveGLS}, 
    \code{objectiveFIML}, \code{msemObjectiveML}, or \code{msemObjectiveGLS}.}
  \item{robust}{if \code{TRUE}, compute robust standard errors or test.}
  \item{fit.indices}{a character vector of ``fit indices'' to report; the allowable values are those given in \bold{Usage}
    above, and vary by the objective function. If the argument isn't given then the fit indices reported are taken
    from the R \code{fit.indices} option; if this option isn't set, then only the AIC and BIC are reported for models fit
    with \code{objectiveML}, \code{objectiveFIML}, or \code{msemObjectiveML}, and no fit indices are reported for
    models fit with \code{objectiveGLS} or \code{msemObjectiveGLS}.}
  \item{k, \ldots}{ignored.}
  \item{digits}{digits to be printed.}
  \item{conf.level}{level for confidence interval for the RMSEA index (default is .9).}
  \item{analytic.se}{use analytic (as opposed to numeric) coefficient standard errors; default is \code{TRUE}
    where analytic standard errors are available if there are no more than
    100 parameters in the model and \code{FALSE} otherwise.}
  \item{saturated}{if \code{TRUE} (the default is \code{FALSE}); compute the log-likelihood (and statistics that
    depend on it) for the saturated model when the objective function is FIML in the presence of missing data. 
    This can be computationally costly.}
  \item{intercept}{the name of the intercept regressor in the raw data, to be used in calculating the
    saturated log-likelihood for the FIML estimator; the default is \code{"Intercept"}.}
  \item{saturated.logLik}{the log-likelihood for the saturated model, as returned by \code{logLik}
    with \code{saturated=TRUE}; if absent, this will be computed and the computation can be time-consuming.}
  \item{iterlim}{iteration limit used by the \code{nlm} optimizer to compute the saturated log-likelihood for
    the FIML estimator with missing data; defaults to \code{1000}.}
}

\author{John Fox \email{jfox@mcmaster.ca} and Jarrett Byrnes}

\references{
See \code{\link{sem}}.
}

\seealso{
\code{\link{sem}}, \code{\link{objective.functions}}, \code{\link{modIndices.objectiveML}} 
}

\keyword{models}
