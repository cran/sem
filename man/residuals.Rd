\name{residuals.sem}
\alias{residuals.sem}
\alias{standardized.residuals}
\alias{standardized.residuals.sem}
\alias{normalized.residuals}
\alias{normalized.residuals.sem}

\title{Residual Covariances for a Structural Equation Model}
\description{
  These functions compute residual covariances, variance-standardized
  residual covariances, and normalized residual covariances
  for the observed variables in a structural-equation
  model fit by \code{sem}.
}
\usage{
\method{residuals}{sem}(object, ...)

standardized.residuals(object, ...)

\method{standardized.residuals}{sem}(object, ...)

normalized.residuals(object, ...)

\method{normalized.residuals}{sem}(object, ...)
}

\arguments{
  \item{object}{an object of class \code{sem} returned by the \code{\link{sem}} function.}
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
    (s_{ij} - c_{ij})/[(c_{ii}c_{jj} - c_{ij}²)/N]^{1/2}}
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
    \dontrun{
# ------------- assumes that Duncan, Haller and Portes peer-influences model
# -------------     has been fit and is in sem.dhp.1

residuals(sem.dhp.1)

##                  ROccAsp        REdAsp       FOccAsp        FEdAsp
##    ROccAsp  4.156103e-07 -2.896368e-07 -2.866110e-02  9.102874e-02
##    REdAsp  -2.896368e-07 -1.541581e-06 -1.094841e-02 -2.379215e-02
##    FOccAsp -2.866110e-02 -1.094841e-02 -3.740356e-06 -9.564103e-07
##    
##    . . .
##    
##                      FIQ       FParAsp
##    ROccAsp  3.467853e-02 -3.187309e-02
##    REdAsp   4.878970e-03 -4.443480e-02
##    FOccAsp -7.354686e-03  2.488120e-02
##    FEdAsp   1.124604e-02 -3.690078e-02
##    RParAsp  2.775558e-17  5.551115e-17
##    RIQ      2.220446e-16  6.938894e-17
##    RSES     0.000000e+00 -1.387779e-17
##    FSES     1.110223e-16 -2.775558e-17
##    FIQ      4.440892e-16  1.110223e-16
##    FParAsp  1.110223e-16  4.440892e-16

normalized.residuals(sem.dhp.1) 

##                  ROccAsp        REdAsp       FOccAsp        FEdAsp
##    ROccAsp  5.330519e-06 -4.455587e-06 -4.898232e-01  1.567678e+00
##    REdAsp  -4.455587e-06 -1.977191e-05 -1.857670e-01 -4.071582e-01
##    FOccAsp -4.898232e-01 -1.857670e-01 -4.797271e-05 -1.460881e-05
##    
##    . . .
##    
##                      FIQ       FParAsp
##    ROccAsp  6.080514e-01 -5.747909e-01
##    REdAsp   8.518738e-02 -8.007295e-01
##    FOccAsp -1.180429e-01  4.374639e-01
##    FEdAsp   1.832159e-01 -6.514685e-01
##    RParAsp  5.019082e-16  1.000322e-15
##    RIQ      3.818356e-15  1.252092e-15
##    RSES     0.000000e+00 -2.506364e-16
##    FSES     1.931472e-15 -5.029583e-16
##    FIQ      5.695780e-15  1.971289e-15
##    FParAsp  1.971289e-15  5.695780e-15
      
standardized.residuals(sem.dhp.1)

##                  ROccAsp        REdAsp       FOccAsp        FEdAsp
##    ROccAsp  4.156103e-07 -2.896368e-07 -2.866110e-02  9.102874e-02
##    REdAsp  -2.896368e-07 -1.541581e-06 -1.094841e-02 -2.379215e-02
##    FOccAsp -2.866110e-02 -1.094841e-02 -3.740356e-06 -9.564103e-07
##    
##    . . .
##    
##                      FIQ       FParAsp
##    ROccAsp  3.467853e-02 -3.187309e-02
##    REdAsp   4.878970e-03 -4.443480e-02
##    FOccAsp -7.354686e-03  2.488120e-02
##    FEdAsp   1.124604e-02 -3.690078e-02
##    RParAsp  2.775558e-17  5.551115e-17
##    RIQ      2.220446e-16  6.938894e-17
##    RSES     0.000000e+00 -1.387779e-17
##    FSES     1.110223e-16 -2.775558e-17
##    FIQ      4.440892e-16  1.110223e-16
##    FParAsp  1.110223e-16  4.440892e-16
    }
}

\keyword{models}
