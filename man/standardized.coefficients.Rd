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
# -------------     has been fit and is in sem.dhp.1     
\dontrun{
standardized.coefficients(sem.dhp.1) 

##           Std. Estimate                     
## 1  gam11  0.210278      RGenAsp <--- RParAsp
## 2  gam12  0.325612      RGenAsp <--- RIQ    
## 3  gam13  0.284855      RGenAsp <--- RSES   
## 4  gam14  0.093702      RGenAsp <--- FSES   
## 5  gam23  0.074576      FGenAsp <--- RSES   
## 6  gam24  0.275763      FGenAsp <--- FSES   
## 7  gam25  0.420558      FGenAsp <--- FIQ    
## 8  gam26  0.192224      FGenAsp <--- FParAsp
## 9  beta12 0.199418      RGenAsp <--- FGenAsp
## 10 beta21 0.217521      FGenAsp <--- RGenAsp
## 11        0.766717      ROccAsp <--- RGenAsp
## 12 lam21  0.814771      REdAsp <--- RGenAsp 
## 13        0.829943      FOccAsp <--- FGenAsp
## 14 lam42  0.771619      FEdAsp <--- FGenAsp 
    }
}

\keyword{models}
