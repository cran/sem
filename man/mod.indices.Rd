\name{mod.indices}
\alias{mod.indices}
\alias{mod.indices.sem}
\alias{print.sem.modind}
\alias{summary.sem.modind}

\title{Modification Indices for Structural Equation Models}
\description{
  \code{mod.indices} calculates modification indices (score tests) 
  and estimated parameter changes for the fixed and constrained 
  parameters in a structural equation model. [Warning: This is a
  preliminary version: the computed values are a bit off the correct
  values.]
}

\usage{
mod.indices(model, ...)

\method{mod.indices}{sem}(model, ...)

\method{print}{sem.modind}(x, n.largest=5, ...)

\method{summary}{sem.modind}(object, round=2, 
    print.matrices=c("both", "par.change", "mod.indices"), ...)
}
    
\arguments{
  \item{model}{an object of class \code{sem}, produced by the \code{\link{sem}} function.}
  \item{object, x}{an object of class \code{sem.modind}, produced by the 
    \code{mod.indices} function.}
  \item{n.largest}{number of modification indices to print in each of the eqn{A} and
    eqn{P} matrices of the RAM model.}
  \item{round}{number of places to the right of the decimal point in printing 
    modification indices.}
  \item{print.matrices}{which matrices to print: estimated changes in the fixed parameters,
    modification indices, or both (the default).}
  \item{\dots}{arguments to be passed down.}
}

\details{
  Modification indices are one-df chi-square score statistics for the fixed and constrained
  parameters in a structural equation model. They may be regarded as an estimate of the improvement
  in the likelihood-ratio chi-square statistic for the model if the corresponding parameter is
  respecified as a free parameter. The \code{mod.indices} function also estimates the change in the
  value of a fixed or constrained parameter if the parameter is respecified as free. When several
  parameters are set equal, modification indices and estimated changes are given for all but the first.
  Modification indices and estimated parameter changes for currently free parameters are given as
  \code{NA}.
  
  The method employed is described in Saris, Satorra, and S\"{o}rbom (1987) and S\"{o}rbom (1989).
}

\value{
  \code{mod.indices} returns an object of class \code{sem.modind} with the following elements:
  \item{mod.A}{modification indices for the elements of the eqn{A} matrix.}
  \item{mod.P}{modification indices for the elements of the eqn{P} matrix.}
  \item{par.A}{estimated parameter changes for the elements of the eqn(A) matrix.}
  \item{par.P}{estimated parameter changes for the elements of the eqn(P) matrix.}
}

\references{
    Sarris, W. E., Satorra, A., and S\"{o}rbom, D. (1987)
    The detection and correction of specification errors in structural
    equation models. Pp. 105--129 in Clogg, C. C. (ed.), \emph{Sociological Methodology 1987.}
    American Sociological Association.
    
    S\"{o}rbom, D. (1989)
    Model modification.
    \emph{Psychometrika} \bold{54}, 371--384.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{
# This example is adapted from the SAS manual

S.wh <- matrix(c(
   11.834,     0,        0,        0,       0,        0,
    6.947,    9.364,     0,        0,       0,        0,
    6.819,    5.091,   12.532,     0,       0,        0,
    4.783,    5.028,    7.495,    9.986,    0,        0,
   -3.839,   -3.889,   -3.841,   -3.625,   9.610,     0,
  -21.899,  -18.831,  -21.748,  -18.775,  35.522,  450.288), 
  6, 6)

model.wh <- matrix(c(
                    'Alienation67   ->  Anomia67',      NA,     1,
                    'Alienation67   ->  Powerless67',   NA,     0.833,
                    'Alienation71   ->  Anomia71',      NA,     1,
                    'Alienation71   ->  Powerless71',   NA,     0.833, 
                    'SES            ->  Education',     NA,     1,     
                    'SES            ->  SEI',           'lamb', NA,
                    'SES            ->  Alienation67',  'gam1', NA,
                    'Alienation67   ->  Alienation71',  'beta', NA,
                    'SES            ->  Alienation71',  'gam2', NA,
                    'Anomia67       <-> Anomia67',      'the1', NA,
                    'Anomia71       <-> Anomia71',      'the1', NA,
                    'Powerless67    <-> Powerless67',   'the2', NA,
                    'Powerless71    <-> Powerless71',   'the2', NA,
                    'Education      <-> Education',     'the3', NA,
                    'SEI            <-> SEI',           'the4', NA,
                    'Anomia67       <-> Anomia71',      'the5', NA,
                    'Powerless67    <-> Powerless71',   'the5', NA,
                    'Alienation67   <-> Alienation67',  'psi1', NA,
                    'Alienation71   <-> Alienation71',  'psi2', NA,
                    'SES            <-> SES',           'phi',  NA), 
                    ncol=3, byrow=TRUE)
                                            
obs.vars.wh <- c('Anomia67','Powerless67','Anomia71','Powerless71','Education','SEI')
rownames(S.wh) <- colnames(S.wh) <- obs.vars.wh

sem.wh <- sem(model.wh, S.wh, 932)

mod.indices(sem.wh)

##   5 largest modification indices, A matrix:
##  Powerless67:Education    Anomia67:Education 
##                 4.8736                3.8027 
##        Powerless67:SES Education:Powerless67 
##                 2.7608                2.4619 
##         Anomia67:SES 
##               2.3122 
##
##    5 largest modification indices, P matrix:
##  Education:Powerless67    Education:Anomia67 
##                 6.4028                4.5398 
##        SES:Powerless67          SES:Anomia67 
##                 2.7608                2.3122 
##        SEI:Powerless67 
##                 1.3185 
}

\keyword{models}
