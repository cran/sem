\name{standardizedCoefficients}
\alias{standardizedCoefficients}
\alias{stdCoef}

\title{Standardized Coefficients for Structural Equation Models}
\description{
  These functions calculate standardized regression coefficients
  for structural equation models. The function \code{stdCoef} is
  simply an abbreviation for \code{standardizedCoefficients}.
}

\usage{
standardizedCoefficients(object, digits = 5, oneheaded=TRUE, twoheaded=TRUE)

stdCoef(...)
}

\arguments{
  \item{object}{an object of class \code{sem} returned by the 
    \code{\link{sem}} function.}
  \item{digits}{number of digits for printed output.}
  \item{oneheaded}{standardize path coefficients? Default is \code{TRUE}.}
  \item{twoheaded}{standardize variances and covariances? Default is \code{TRUE}.}
  \item{...}{arguments to pass to \code{standardized.coefficients}.}
}

\value{
  Returns a data frame with the coefficients, labelled 
  both by parameter names and by arrows in the path diagram for the
  model.
}

\references{
    Bollen, K. A. (1989) 
    \emph{Structural Equations With Latent Variables.} Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca} and Adam Kramer}

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
standardizedCoefficients(sem.dhp) 
	}
# The following example can be executed via example():

etc <- file.path(.path.package(package="sem")[1], "etc") # path to data and model files
   
(R.DHP <- readMoments(file=file.path(etc, "R-DHP.txt"),
				diag=FALSE, names=c("ROccAsp", "REdAsp", "FOccAsp", 
                "FEdAsp", "RParAsp", "RIQ", "RSES", "FSES", "FIQ", "FParAsp")))
(model.dhp <- specifyModel(file=file.path(etc, "model-DHP.txt")))
(sem.dhp <- sem(model.dhp, R.DHP, 329,
    fixed.x=c('RParAsp', 'RIQ', 'RSES', 'FSES', 'FIQ', 'FParAsp')))
standardizedCoefficients(sem.dhp) 
}

\keyword{models}
