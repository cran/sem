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
residuals.sem(object)

standardized.residuals(object, ...)

standardized.residuals.sem(object)

normalized.residuals(object, ...)

normalized.residuals.sem(object)
}

\arguments{
  \item{object}{an object of class \code{sem} returned by the \code{\link{sem}} function.}
}

\details{
  Residuals are defined as \eqn{S - C}, where \eqn{S} is the sample covariance matrix
  of the observed variables and \eqn{C} is the model-reproduced covariance matrix.
  The \emph{standardized} residual covariance for a pair of variables divides the
  residual covariance by the product of the sample standard deviations of the 
  two variables, \eqn{(s_{ij} - c_{ij})/(s_{ii}s_{jj})^{1/2}}. The \emph{normalized} residual
  is given by
  \deqn{\frac{s_{ij}-c_{ij}} {[(c_{ii}c_{jj}-c_{ij}^2)/N]^{1/2}}}{%
    (s_{ij} - c_{ij})/[(c_{ii}c_{jj} - c_{ij}²)/N]^{1/2}}
  where \eqn{N} is the number of observations.
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
# -------------     has been fit and is in sem.dhp

residuals(sem.dhp)

##                ROccAsp        REdAsp       FOccAsp        FEdAsp
##    ROccAsp -2.012270e-07  1.665266e-07 -2.866299e-02  9.102173e-02
##    REdAsp   1.665266e-07  2.032838e-06 -1.094733e-02 -2.379674e-02
##    FOccAsp -2.866299e-02 -1.094733e-02  3.517460e-06  1.093859e-06
##    
##    . . .
##    
##                    FIQ       FParAsp
##    ROccAsp  3.467268e-02 -3.187624e-02
##    REdAsp   4.875055e-03 -4.443721e-02
##    FOccAsp -7.350068e-03  2.488664e-02
##    FEdAsp   1.124255e-02 -3.689947e-02
##    RParAsp  1.387779e-17  5.551115e-17
##    RIQ     -5.551115e-17  0.000000e+00
##    RSES     5.551115e-17  1.387779e-17
##    FSES    -5.551115e-17  0.000000e+00
##    FIQ      0.000000e+00  2.775558e-17
##    FParAsp  2.775558e-17  4.440892e-16

normalized.residuals(sem.dhp) 

##                ROccAsp        REdAsp       FOccAsp        FEdAsp
##    ROccAsp -2.580888e-06  2.561743e-06 -4.898566e-01  1.567552e+00
##    REdAsp   2.561743e-06  2.607273e-05 -1.857497e-01 -4.072364e-01
##    FOccAsp -4.898566e-01 -1.857497e-01  4.511424e-05  1.670834e-05
##    
##    . . .
##                    FIQ       FParAsp
##    ROccAsp  6.079479e-01 -5.748474e-01
##    REdAsp   8.511907e-02 -8.007742e-01
##    FOccAsp -1.179693e-01  4.375617e-01
##    FEdAsp   1.831587e-01 -6.514448e-01
##    RParAsp  2.509541e-16  1.000322e-15
##    RIQ     -9.545890e-16  0.000000e+00
##    RSES     9.812183e-16  2.506364e-16
##    FSES    -9.657360e-16  0.000000e+00
##    FIQ      0.000000e+00  4.928223e-16
##    FParAsp  4.928223e-16  5.695780e-15

standardized.residuals(sem.dhp)

##                ROccAsp        REdAsp       FOccAsp        FEdAsp
##    ROccAsp -2.012270e-07  1.665266e-07 -2.866299e-02  9.102173e-02
##    REdAsp   1.665266e-07  2.032838e-06 -1.094733e-02 -2.379674e-02
##    FOccAsp -2.866299e-02 -1.094733e-02  3.517460e-06  1.093859e-06
##    
##    . . .
##    
##                    FIQ       FParAsp
##    ROccAsp  3.467268e-02 -3.187624e-02
##    REdAsp   4.875055e-03 -4.443721e-02
##    FOccAsp -7.350068e-03  2.488664e-02
##    FEdAsp   1.124255e-02 -3.689947e-02
##    RParAsp  1.387779e-17  5.551115e-17
##    RIQ     -5.551115e-17  0.000000e+00
##    RSES     5.551115e-17  1.387779e-17
##    FSES    -5.551115e-17  0.000000e+00
##    FIQ      0.000000e+00  2.775558e-17
##    FParAsp  2.775558e-17  4.440892e-16
    }
}

\keyword{models}
