\name{ram}
\alias{ram}

\title{RAM Matrix for a Structural-Equation Model}
\description{
  Print the labelled RAM definition matrix for a structural-equation
  model fit by \code{sem}.
}

\usage{
ram(object, digits=5, startvalues=F)
}


\arguments{
  \item{object}{an object of class \code{sem} returned by the \code{sem} function.}
  \item{digits}{number of digits for printed output.}
  \item{startvalues}{if \code{TRUE}, start values for parameters are printed;
    otherwise, the parameter estimates are printed; the default is \code{FALSE}.}
}


}
\value{
  A data frame containing the labelled RAM definition matrix, which is normally
  just printed.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{

# ------------- assumes that Duncan, Haller and Portes peer-influences model
# -------------     has been fit and is in sem.dhp

    \dontrun{
ram(sem.dhp)

##        heads to from parameter  estimate                arrow
##            1  1   11         0  1.000000 ROccAsp <--- RGenAsp
## lam21      1  2   11         1  1.062673  REdAsp <--- RGenAsp
##            1  3   12         0  1.000000 FOccAsp <--- FGenAsp
## lam42      1  4   12         2  0.929732  FEdAsp <--- FGenAsp
## gam11      1 11    5         3  0.161220 RGenAsp <--- RParAsp
## gam12      1 11    6         4  0.249647     RGenAsp <--- RIQ
## gam13      1 11    7         5  0.218402    RGenAsp <--- RSES
## gam14      1 11    8         6  0.071836    RGenAsp <--- FSES
## gam23      1 12    7         7  0.061879    FGenAsp <--- RSES
## gam24      1 12    8         8  0.228863    FGenAsp <--- FSES
## gam25      1 12    9         9  0.349030     FGenAsp <--- FIQ
## gam26      1 12   10        10  0.159529 FGenAsp <--- FParAsp
## beta12     1 11   12        11  0.184245 RGenAsp <--- FGenAsp
## beta21     1 12   11        12  0.235502 FGenAsp <--- RGenAsp
## theta1     2  1    1        13  0.412143 ROccAsp <--> ROccAsp
## theta2     2  2    2        14  0.336146   REdAsp <--> REdAsp
## theta3     2  3    3        15  0.311197 FOccAsp <--> FOccAsp
## theta4     2  4    4        16  0.404601   FEdAsp <--> FEdAsp
## psi11      2 11   11        17  0.280987 RGenAsp <--> RGenAsp
## psi22      2 12   12        18  0.263832 FGenAsp <--> FGenAsp
## psi12      2 11   12        19 -0.022620 RGenAsp <--> FGenAsp
##            2  5    5         0  1.000000 RParAsp <--> RParAsp
##            2  6    5         0  0.183900     RIQ <--> RParAsp
##            2  6    6         0  1.000000         RIQ <--> RIQ
##            2  7    5         0  0.048900    RSES <--> RParAsp
##            2  7    6         0  0.222000        RSES <--> RIQ
##            2  7    7         0  1.000000       RSES <--> RSES
##            2  8    5         0  0.018600    FSES <--> RParAsp
##            2  8    6         0  0.186100        FSES <--> RIQ
##            2  8    7         0  0.270700       FSES <--> RSES
##            2  8    8         0  1.000000       FSES <--> FSES
##            2  9    5         0  0.078200     FIQ <--> RParAsp
##            2  9    6         0  0.335500         FIQ <--> RIQ
##            2  9    7         0  0.230200        FIQ <--> RSES
##            2  9    8         0  0.295000        FIQ <--> FSES
##            2  9    9         0  1.000000         FIQ <--> FIQ
##            2 10    5         0  0.114700 FParAsp <--> RParAsp
##            2 10    6         0  0.102100     FParAsp <--> RIQ
##            2 10    7         0  0.093100    FParAsp <--> RSES
##            2 10    8         0 -0.043800    FParAsp <--> FSES
##            2 10    9         0  0.208700     FParAsp <--> FIQ
##            2 10   10         0  1.000000 FParAsp <--> FParAsp
    }
}
\keyword{models}
