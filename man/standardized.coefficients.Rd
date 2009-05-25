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
standardized.coefficients(object, digits = 5, oneheaded=TRUE, twoheaded=TRUE)

std.coef(...)
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
  Returns a data frame with the standardized coefficients, labelled 
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

# ------------- assumes that Duncan, Haller and Portes peer-influences model
# -------------     has been fit and is in sem.dhp.1     
\dontrun{
standardized.coefficients(sem.dhp.1) 

##           Std. Estimate                     
          Std. Estimate                     
## 1  gam11   0.210278     RGenAsp <--- RParAsp
## 2  gam12   0.325612     RGenAsp <--- RIQ    
## 3  gam13   0.284855     RGenAsp <--- RSES   
## 4  gam14   0.093702     RGenAsp <--- FSES   
## 5  gam23   0.074576     FGenAsp <--- RSES   
## 6  gam24   0.275763     FGenAsp <--- FSES   
## 7  gam25   0.420558     FGenAsp <--- FIQ    
## 8  gam26   0.192224     FGenAsp <--- FParAsp
## 9  beta12  0.199418     RGenAsp <--- FGenAsp
## 10 beta21  0.217521     FGenAsp <--- RGenAsp
## 11         0.766717     ROccAsp <--- RGenAsp
## 12 lam21   0.814771     REdAsp <--- RGenAsp 
## 13         0.829943     FOccAsp <--- FGenAsp
## 14 lam42   0.771619     FEdAsp <--- FGenAsp 
## 15 ps11    0.477987     RGenAsp <--> RGenAsp
## 16 ps22    0.383036     FGenAsp <--> FGenAsp
## 17 ps12   -0.035518     FGenAsp <--> RGenAsp
## 18 theta1  0.412144     ROccAsp <--> ROccAsp
## 19 theta2  0.336148     REdAsp <--> REdAsp  
## 20 theta3  0.311194     FOccAsp <--> FOccAsp
## 21 theta4  0.404603     FEdAsp <--> FEdAsp  
## 22         1.000000     RParAsp <--> RParAsp
## 23         0.183900     RIQ <--> RParAsp    
## 24         1.000000     RIQ <--> RIQ        
## 25         0.048900     RSES <--> RParAsp   
## 26         0.222000     RSES <--> RIQ       
## 27         1.000000     RSES <--> RSES      
## 28         0.018600     FSES <--> RParAsp   
## 29         0.186100     FSES <--> RIQ       
## 30         0.270700     FSES <--> RSES      
## 31         1.000000     FSES <--> FSES      
## 32         0.078200     FIQ <--> RParAsp    
## 33         0.335500     FIQ <--> RIQ        
## 34         0.230200     FIQ <--> RSES       
## 35         0.295000     FIQ <--> FSES       
## 36         1.000000     FIQ <--> FIQ        
## 37         0.114700     FParAsp <--> RParAsp
## 38         0.102100     FParAsp <--> RIQ    
## 39         0.093100     FParAsp <--> RSES   
## 40        -0.043800     FParAsp <--> FSES   
## 41         0.208700     FParAsp <--> FIQ    
## 42         1.000000     FParAsp <--> FParAsp
    }
}

\keyword{models}
