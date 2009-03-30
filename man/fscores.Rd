\name{fscores}
\alias{fscores}
\alias{fscores.sem}

\title{Factor Scores for Latent Variables}
\description{
  Calculate factor scores or factor-score coefficients for the latent variables
  in a structural-equation model.
}
\usage{
fscores(model, ...)

\method{fscores}{sem}(model, data, center = TRUE, scale = FALSE, ...)
}

\arguments{
  \item{model}{an object of class \code{"sem"}, produced by the \code{\link{sem}}
    function.}
  \item{data}{an optional numeric data frame or matrix containing the observed variables
    in the model; if present, the estimated factor scores are returned; if absent, the
    factor-score \emph{coefficients} are returned.}
  \item{center}{if \code{TRUE}, the default, the means of the observed variables are
    subtracted prior to computing factor scores. One would normally use this option
    if the model is estimated from a covariance or correlation matrix among the
    observed variables.}
  \item{scale}{if \code{TRUE}, the possibly centered variables are divided by their
    room-mean-squares; the default is \code{FALSE}.
    One would normally use this option if the model is estimated
    from a correlation matrix among the observed variables. Centering and scaling 
    are performed by the \code{\link[base]{scale}} function.}
  \item{\dots}{arguments to pass down.}
}
\details{
  Factor-score coefficients are computed by the \dQuote{regression} method as
  \eqn{B = C^{-1} C^{*}}{C^-1 C*}, where \eqn{C} is the model-implied covariance or
  moment matrix among the observed variables and  \eqn{C^{*}}{C*} is the matrix
  of model-implied covariances or moments between the observed and latent variables.
}
\value{
  Either a matrix of estimated factor scores (if the \code{data} argument is
  supplied) or a matrix of factor-score coefficients (otherwise).
}

\references{
    Bollen, K. A. (1989) 
    \emph{Structural Equations With Latent Variables.} Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}, \code{\link[base]{scale}}}

\examples{
	   \dontrun{

S.wh <- read.moments(names=c('Anomia67','Powerless67','Anomia71',
                                    'Powerless71','Education','SEI'))
   11.834                                    
    6.947    9.364                            
    6.819    5.091   12.532                    
    4.783    5.028    7.495    9.986            
   -3.839   -3.889   -3.841   -3.625   9.610     
  -21.899  -18.831  -21.748  -18.775  35.522  450.288

# This model in the SAS manual for PROC CALIS

model.wh.1 <- specify.model()
    Alienation67   ->  Anomia67,      NA,     1
    Alienation67   ->  Powerless67,   NA,     0.833
    Alienation71   ->  Anomia71,      NA,     1
    Alienation71   ->  Powerless71,   NA,     0.833 
    SES            ->  Education,     NA,     1     
    SES            ->  SEI,           lamb,   NA
    SES            ->  Alienation67,  gam1,   NA
    Alienation67   ->  Alienation71,  beta,   NA
    SES            ->  Alienation71,  gam2,   NA
    Anomia67       <-> Anomia67,      the1,   NA
    Anomia71       <-> Anomia71,      the1,   NA
    Powerless67    <-> Powerless67,   the2,   NA
    Powerless71    <-> Powerless71,   the2,   NA
    Education      <-> Education,     the3,   NA
    SEI            <-> SEI,           the4,   NA
    Anomia67       <-> Anomia71,      the5,   NA
    Powerless67    <-> Powerless71,   the5,   NA
    Alienation67   <-> Alienation67,  psi1,   NA
    Alienation71   <-> Alienation71,  psi2,   NA
    SES            <-> SES,           phi,    NA
    
                        
sem.wh.1 <- sem(model.wh.1, S.wh, 932)

fscores(sem.wh.1)

##                Alienation67 Alienation71         SES
##    Anomia67     0.413112363  0.048268330 -0.05212632
##    Powerless67  0.345402079  0.040014780 -0.04355578
##    Anomia71     0.052663484  0.430618716 -0.03999218
##    Powerless71  0.043704122  0.360044434 -0.03339943
##    Education   -0.074921670 -0.063969383  0.50571037
##    SEI         -0.004638977 -0.003960837  0.03131242

   }
}

\keyword{models}
\keyword{manip}
