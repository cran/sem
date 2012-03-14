\name{modIndices}
\alias{modIndices}
\alias{modIndices.objectiveML}
\alias{modIndices.msemObjectiveML}
\alias{print.modIndices}
\alias{print.msemModIndices}
\alias{summary.modIndices}
\alias{summary.msemModIndices}

\title{Modification Indices for Structural Equation Models}
\description{
  \code{mod.indices} calculates modification indices (score tests) 
  and estimated parameter changes for the fixed and constrained 
  parameters in a structural equation model fit by multinormal maximum likelihood.
}

\usage{
\method{modIndices}{objectiveML}(model, duplicated, deviance=NULL, ...)
\method{modIndices}{msemObjectiveML}(model, ...)

\method{print}{modIndices}(x, n.largest=5, ...)
\method{print}{msemModIndices}(x, ...)

\method{summary}{modIndices}(object, round=2, 
    print.matrices=c("both", "par.change", "mod.indices"), ...)
\method{summary}{msemModIndices}(object, ...)
}
    
\arguments{
  \item{model}{an object of class \code{objectiveML} or \code{msemObjectiveML},
    produced by the \code{\link{sem}} function.}
  \item{object, x}{an object of class \code{modIndices} or \code{msemModIndices}, produced by the 
    \code{modIndices} function.}
  \item{n.largest}{number of modification indices to print in each of the \eqn{A} and
    \eqn{P} matrices of the RAM model.}
  \item{round}{number of places to the right of the decimal point in printing 
    modification indices.}
  \item{print.matrices}{which matrices to print: estimated changes in the fixed parameters,
    modification indices, or both (the default).}
  \item{duplicated, deviance}{for internal use.}
  \item{\dots}{arguments to be passed down.}
}

\details{
  Modification indices are one-df chi-square score (\dQuote{Lagrange-multiplier}) test statistics for the fixed and constrained
  parameters in a structural equation model. They may be regarded as an estimate of the improvement
  in the likelihood-ratio chi-square statistic for the model if the corresponding parameter is
  respecified as a free parameter. The \code{modIndices} function also estimates the change in the
  value of a fixed or constrained parameter if the parameter is respecified as free. When several
  parameters are set equal, modification indices and estimated changes are given for all but the first.
  Modification indices and estimated parameter changes for currently free parameters are given as
  \code{NA}.
  
  The method employed is described in Saris, Satorra, and Sorbom (1987) and Sorbom (1989).
}

\value{
  \code{modIndices} returns an object of class \code{modIndices} with the following elements:
  \item{mod.A}{modification indices for the elements of the \eqn{A} matrix.}
  \item{mod.P}{modification indices for the elements of the \eqn{P} matrix.}
  \item{par.A}{estimated parameter changes for the elements of the \eqn{A} matrix.}
  \item{par.P}{estimated parameter changes for the elements of the \eqn{P} matrix.}
}

\references{
    Sarris, W. E., Satorra, A., and Sorbom, D. (1987)
    The detection and correction of specification errors in structural
    equation models. Pp. 105--129 in Clogg, C. C. (ed.), \emph{Sociological Methodology 1987.}
    American Sociological Association.
    
    Sorbom, D. (1989)
    Model modification.
    \emph{Psychometrika} \bold{54}, 371--384.
}

\author{John Fox \email{jfox@mcmaster.ca} and Michael Culbertson}

\seealso{\code{\link{sem}}}

\examples{
# In the first example, readMoments() and specifyModel() read from the
# input stream. This example cannot be executed via example() but can be entered
# at the command prompt. The example is repeated using file input;
# this example can be executed via example(). 
	\dontrun{
# This example is adapted from the SAS manual

S.wh <- readMoments(names=c('Anomia67','Powerless67','Anomia71',
                                    'Powerless71','Education','SEI'))
   11.834                                    
    6.947    9.364                            
    6.819    5.091   12.532                    
    4.783    5.028    7.495    9.986            
   -3.839   -3.889   -3.841   -3.625   9.610     
  -21.899  -18.831  -21.748  -18.775  35.522  450.288

model.wh <- specifyModel()
    Alienation67   ->  Anomia67,      NA,   1
    Alienation67   ->  Powerless67,   NA,   0.833
    Alienation71   ->  Anomia71,      NA,   1
    Alienation71   ->  Powerless71,   NA,   0.833
    SES            ->  Education,     NA,   1     
    SES            ->  SEI,           lamb, NA
    SES            ->  Alienation67,  gam1, NA
    Alienation67   ->  Alienation71,  beta, NA
    SES            ->  Alienation71,  gam2, NA
    Anomia67       <-> Anomia67,      the1, NA
    Anomia71       <-> Anomia71,      the1, NA
    Powerless67    <-> Powerless67,   the2, NA
    Powerless71    <-> Powerless71,   the2, NA
    Education      <-> Education,     the3, NA
    SEI            <-> SEI,           the4, NA
    Anomia67       <-> Anomia71,      the5, NA
    Powerless67    <-> Powerless71,   the5, NA
    Alienation67   <-> Alienation67,  psi1, NA
    Alienation71   <-> Alienation71,  psi2, NA
    SES            <-> SES,           phi,  NA

sem.wh <- sem(model.wh, S.wh, 932)
modIndices(sem.wh)
	}
	
# The following example can be executed via example():

etc <- file.path(.path.package(package="sem")[1], "etc") # path to data and model files

(S.wh <- readMoments(file=file.path(etc, "S-Wheaton.txt"),
					names=c('Anomia67','Powerless67','Anomia71',
                            'Powerless71','Education','SEI')))
(model.wh <- specifyModel(file=file.path(etc, "model-Wheaton-1.txt")))                    
(sem.wh <- sem(model.wh, S.wh, 932))
modIndices(sem.wh)
}

\keyword{models}
