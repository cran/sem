\name{optimizers}

\alias{optimizers}
\alias{optimizerSem}
\alias{optimizerNlm}
\alias{optimizerOptim}
\alias{optimizerNlminb}
\alias{optimizerMsem}
\alias{msemOptimizerNlm}

\title{sem Optimizers}

\description{
The default optimizer used by \code{\link{sem}} is \code{optimizerSem}, which employs compiled code and is integrated with
the \code{\link{objectiveML}} and \code{\link{objectiveGLS}} objective functions; \code{optimizerSem} is a modified
version of the standard R \code{nlm} optimizer.
The other functions call optimizers (\code{\link[stats]{nlm}}, \code{\link[stats]{optim}}, or \code{\link[stats]{nlminb}}), 
to fit structural equation models, and are called by the \code{\link{sem}} function. 
The user would not normally call these functions directly, but rather supply one of them in the \code{optimizer} argument to
\code{sem}. Users may also write them own optimizer functions. \code{msemOptimizerNlm} is for fitting multigroup models.
}

\usage{
optimizerSem(start, objective=objectiveML,  
	gradient=TRUE, maxiter, debug, par.size, model.description, warn, ...)
	
optimizerMsem(start, objective=msemObjectiveML, gradient=TRUE,
	maxiter, debug, par.size, model.description, warn=FALSE, ...)
	
optimizerNlm(start, objective=objectiveML, gradient=TRUE, 
	maxiter, debug, par.size, model.description, warn, ...)
	
optimizerOptim(start, objective=objectiveML, gradient=TRUE, 
	maxiter, debug, par.size, model.description, warn, method="CG", ...)

optimizerNlminb(start, objective=objectiveML, gradient=TRUE, maxiter, 
	debug, par.size, model.description, warn, ...)

msemOptimizerNlm(start, objective=msemObjectiveML, gradient=TRUE,
		maxiter, debug, par.size, model.description, warn=FALSE, ...)
}

\arguments{

  \item{start}{a vector of start values for the parameters.}
  
  \item{objective}{the objective function to be optimized; see \link{objective.functions}.}
  
  \item{gradient}{\code{TRUE} if an analytic gradient is to be used (if one is available).}
  
  \item{maxiter}{the maximum number of iterations allowed.}
  
  \item{debug}{\code{TRUE} to show the iteration history and other available information about the optimization.}
  
  \item{par.size}{\code{"startvalues"} to have the optimizer scale the problem according to the magitudes of the start values (ignored by \code{optimizerNlminb}).}
  
  \item{model.description}{a list with elements describing the structural-equation model (see the code for details).}
  
  \item{warn}{if \code{FALSE}, suppress warnings during the optimization.}
  
  \item{method}{the method to be employed by the \code{\link[stats]{optim}} optimizer; the default is \code{"CG"} (conjugate-gradient).}
  
  \item{\dots}{additional arguments for the \code{nlm}, \code{optim}, or \code{nlminb} optimizer.}

}

\value{
An object of class \code{"semResult"}, with elements:

\item{convergence}{\code{TRUE} if the optimization apparently converged.}

\item{iterations}{the number of iterations required.}

\item{par}{the vector of parameter estimates.}

\item{vcov}{the estimated covariance matrix of the parameter estimates, based on a numeric Hessian; not supplied by \code{optimizerNlminb}.}

\item{criterion}{the optimized value of the objective function.}

\item{C}{the model-implied covariance or moment matrix at the parameter estimates.}

\item{A}{the estimated \eqn{A} matrix.}

\item{P}{the estimated \eqn{P} matrix.}

}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}, \code{\link{objective.functions}}, \code{\link[stats]{nlm}}, \code{\link[stats]{optim}}, \code{\link[stats]{nlminb}}}

\keyword{models}
