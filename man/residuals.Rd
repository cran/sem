\name{residuals.sem}
\alias{residuals.sem}
\alias{residuals.msem}
\alias{standardizedResiduals}
\alias{standardizedResiduals.sem}
\alias{standardizedResiduals.msem}
\alias{normalizedResiduals}
\alias{normalizedResiduals.objectiveML}
\alias{normalizedResiduals.objectiveGLS}
\alias{normalizedResiduals.msemObjectiveML}

\title{Residual Covariances for a Structural Equation Model}
\description{
  These functions compute residual covariances, variance-standardized
  residual covariances, and normalized residual covariances
  for the observed variables in a structural-equation
  model fit by \code{sem}.
}
\usage{
\method{residuals}{sem}(object, ...)
\method{residuals}{msem}(object, ...)

\method{standardizedResiduals}{sem}(object, ...)
\method{standardizedResiduals}{msem}(object, ...)

\method{normalizedResiduals}{objectiveML}(object, ...)
\method{normalizedResiduals}{objectiveGLS}(object, ...)
\method{normalizedResiduals}{msemObjectiveML}(object, ...)
}

\arguments{
  \item{object}{an object inheriting from class \code{sem} or \code{msem} returned by the \code{\link{sem}} function.}
  \item{...}{not for the user.}
}

\details{
  Residuals are defined as \eqn{S - C}, where \eqn{S} is the sample covariance matrix
  of the observed variables and \eqn{C} is the model-reproduced covariance matrix.
  The \emph{standardized} residual covariance for a pair of variables divides the
  residual covariance by the product of the sample standard deviations of the 
  two variables, \eqn{(s_{ij} - c_{ij})/(s_{ii}s_{jj})^{1/2}}. The \emph{normalized} residual
  is given by
  \deqn{\frac{s_{ij}-c_{ij}} {[(c_{ii}c_{jj}-c_{ij}^2)/N^{*}]^{1/2}}}{%
    (s[ij] - c[ij])/[(c[ii]c[ii] + c[ij]^2)/N*]^[1/2]}
  where \eqn{N^{*}} is the number of observations minus one if the model is fit to a
  covariance matrix, or the number of observations if it is fit to a raw moment matrix.
}

\value{
  Each function returns a matrix of residuals.
}

\references{
    Bollen, K. A. (1989) 
    \emph{Structural Equations With Latent Variables.} Wiley.
    }
    
\author{John Fox \email{jfox@mcmaster.ca}}


\seealso{\code{\link{sem}}}

\examples{
# In the first example, readMoments() and specifyModel() read from the
# input stream. This example cannot be executed via example() but can be entered
# at the command prompt. The example is repeated using file input;
# this example can be executed via example(). 
    \dontrun{
# Duncan, Haller, and Portes peer-influences model

R.DHP <- readMoments(diag=FALSE, names=c("ROccAsp", "REdAsp", "FOccAsp", 
                "FEdAsp", "RParAsp", "RIQ", "RSES", "FSES", "FIQ", "FParAsp"))
    .6247     
    .3269  .3669       
    .4216  .3275  .6404
    .2137  .2742  .1124  .0839
    .4105  .4043  .2903  .2598  .1839
    .3240  .4047  .3054  .2786  .0489  .2220
    .2930  .2407  .4105  .3607  .0186  .1861  .2707
    .2995  .2863  .5191  .5007  .0782  .3355  .2302  .2950
    .0760  .0702  .2784  .1988  .1147  .1021  .0931 -.0438  .2087
            
model.dhp <- specifyModel()
    RParAsp  -> RGenAsp, gam11,  NA
    RIQ      -> RGenAsp, gam12,  NA
    RSES     -> RGenAsp, gam13,  NA
    FSES     -> RGenAsp, gam14,  NA
    RSES     -> FGenAsp, gam23,  NA
    FSES     -> FGenAsp, gam24,  NA
    FIQ      -> FGenAsp, gam25,  NA
    FParAsp  -> FGenAsp, gam26,  NA
    FGenAsp  -> RGenAsp, beta12, NA
    RGenAsp  -> FGenAsp, beta21, NA
    RGenAsp  -> ROccAsp,  NA,     1
    RGenAsp  -> REdAsp,  lam21,  NA
    FGenAsp  -> FOccAsp,  NA,     1
    FGenAsp  -> FEdAsp,  lam42,  NA
    RGenAsp <-> RGenAsp, ps11,   NA
    FGenAsp <-> FGenAsp, ps22,   NA
    RGenAsp <-> FGenAsp, ps12,   NA
    ROccAsp <-> ROccAsp, theta1, NA
    REdAsp  <-> REdAsp,  theta2, NA
    FOccAsp <-> FOccAsp, theta3, NA
    FEdAsp  <-> FEdAsp,  theta4, NA

sem.dhp <- sem(model.dhp, R.DHP, 329,
    fixed.x=c('RParAsp', 'RIQ', 'RSES', 'FSES', 'FIQ', 'FParAsp'))
residuals(sem.dhp)
normalizedResiduals(sem.dhp) 
standardizedResiduals(sem.dhp) # same as residuals because model is fit to correlations
    }
# The following example can be executed via example():

etc <- file.path(.path.package(package="sem")[1], "etc") # path to data and model files
   
(R.DHP <- readMoments(file=file.path(etc, "R-DHP.txt"),
				diag=FALSE, names=c("ROccAsp", "REdAsp", "FOccAsp", 
                "FEdAsp", "RParAsp", "RIQ", "RSES", "FSES", "FIQ", "FParAsp")))
(model.dhp <- specifyModel(file=file.path(etc, "model-DHP.txt")))
(sem.dhp <- sem(model.dhp, R.DHP, 329,
    fixed.x=c('RParAsp', 'RIQ', 'RSES', 'FSES', 'FIQ', 'FParAsp')))

residuals(sem.dhp)

normalizedResiduals(sem.dhp) 

standardizedResiduals(sem.dhp)  # same as residuals because model is fit to correlations
}

\keyword{models}
