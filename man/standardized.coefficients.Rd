\name{standardized.coefficients}
\alias{standardized.coefficients}
\alias{std.coef}

\title{Standardized Coefficients for Structural Equation Models}
\description{
  These functions calculate standardized regression coefficients
  for structural equation models. The function \code{std.coef} is
  simply an abbreviation for \code{standardized.coefficients}.
}

\usage{
standardized.coefficients(object, digits=5)

std.coef(...)
}

\arguments{
  \item{object}{an object of class \code{sem} returned by the 
    \code{\link{sem}} function.}
  \item{digits}{number of digits for printed output.}
  \item{...}{arguments to pass to \code{standardized.coefficients}.}
}

\value{
  Returns a data frame with the standardized coefficients, labelled 
  both by parameter names and by arrows in the path diagram for the
  model.
}

\references{
    Bollen, K. A. (1989) 
    \emph{Structural Equations With Latent Variables.} Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{

# ------------- assumes that Duncan, Haller and Portes peer-influences model
# -------------     has been fit and is in sem.dhp
    \dontrun{
standardized.coefficients(sem.dhp) 

##            Std. Estimate                     
##               0.766718 ROccAsp <--- RGenAsp
##    lam21      0.814772  REdAsp <--- RGenAsp
##               0.829941 FOccAsp <--- FGenAsp
##    lam42      0.771621  FEdAsp <--- FGenAsp
##    gam11      0.210273 RGenAsp <--- RParAsp
##    gam12      0.325604     RGenAsp <--- RIQ
##    gam13      0.284853    RGenAsp <--- RSES
##    gam14      0.093693    RGenAsp <--- FSES
##    gam23      0.074559    FGenAsp <--- RSES
##    gam24      0.275758    FGenAsp <--- FSES
##    gam25      0.420549     FGenAsp <--- FIQ
##    gam26      0.192218 FGenAsp <--- FParAsp
##    beta12     0.199437 RGenAsp <--- FGenAsp
##    beta21     0.217562 FGenAsp <--- RGenAsp
    }
}

\keyword{models}
