\name{fscores}
\alias{fscores}
\alias{fscores.sem}
\alias{fscores.msem}

\title{Factor Scores for Latent Variables}
\description{
  Calculate factor scores or factor-score coefficients for the latent variables
  in a structural-equation model.
}
\usage{
\method{fscores}{sem}(model, data=model$data, center=TRUE, scale=FALSE, ...)
\method{fscores}{msem}(model, data=model$data, center=TRUE, scale=FALSE, ...)
}

\arguments{
  \item{model}{an object of class \code{"sem"} or \code{"msem"}, produced by the \code{\link{sem}}
    function.}
  \item{data}{an optional numeric data frame or matrix containing the observed variables
    in the model; if not \code{NULL}, the estimated factor scores are returned; if \code{NULL}, the
    factor-score \emph{coefficients} are returned. The default is the \code{data} element of \code{model},
    which is non-\code{NULL} if the model was fit to a data set rather than a covariance or moment matrix.}
  \item{center}{if \code{TRUE}, the default, the means of the observed variables are
    subtracted prior to computing factor scores. One would normally use this option
    if the model is estimated from a covariance or correlation matrix among the
    observed variables.}
  \item{scale}{if \code{TRUE}, the possibly centered variables are divided by their
    root-mean-squares; the default is \code{FALSE}.
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
  supplied) or a matrix of factor-score coefficients (otherwise). In the case of an \code{"msem"}
  argument, a list of matrices is returned.
}

\references{
    Bollen, K. A. (1989) 
    \emph{Structural Equations With Latent Variables.} Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}, \code{\link[base]{scale}}}

\examples{

# In the first example, readMoments() and specifyModel() read from the
# input stream. This example cannot be executed via example() but can be entered
# at the command prompt. The example is repeated using file input;
# this example can be executed via example(). 
	   \dontrun{

S.wh <- readMoments(names=c('Anomia67','Powerless67','Anomia71',
                                    'Powerless71','Education','SEI'))
   11.834                                    
    6.947    9.364                            
    6.819    5.091   12.532                    
    4.783    5.028    7.495    9.986            
   -3.839   -3.889   -3.841   -3.625   9.610     
  -21.899  -18.831  -21.748  -18.775  35.522  450.288

# This model in the SAS manual for PROC CALIS

model.wh.1 <- specifyModel()
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
   }

# The following example can be executed via example():

etc <- file.path(.path.package(package="sem")[1], "etc") # path to data and model files
   
(S.wh <- readMoments(file=file.path(etc, "S-Wheaton.txt"),
					names=c('Anomia67','Powerless67','Anomia71',
                            'Powerless71','Education','SEI')))
(model.wh.1 <- specifyModel(file=file.path(etc, "model-Wheaton-1.txt")))        
(sem.wh.1 <- sem(model.wh.1, S.wh, 932))
fscores(sem.wh.1)
}

\keyword{models}
\keyword{manip}
