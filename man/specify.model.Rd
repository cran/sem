\name{specify.model}
\alias{specify.model}
\alias{print.mod}
\title{Specify a Structural Equation Model}
\description{
  Create the RAM specification of a structural equation model.
}
\usage{
specify.model(file = "")

\method{print}{mod}(x, ...)
}

\arguments{
  \item{file}{The (quoted) file from which to read the model specification,
    including the path to the file if it is not in the current directory. If
    \code{""} (the default), then the specification is read from the standard
    input stream, and is terminated by a blank line.}
  \item{x}{An object of class \code{mod} to print, as produced by
    \code{specify.model}.}
  \item{...}{Ignored.}
}

\details{
  Each line of the RAM specification consists of three (unquoted) entries, 
  separated by commas:
  
    \describe{
      \item{1. Arrow specification:}{This is a simple formula, of the form
        \code{A -> B} or, equivalently, \code{B <- A} for a regression
        coefficient (i.e., a single-headed or directional arrow); 
        \code{A <-> A} for a variance or \code{A <-> B} for a covariance
        (i.e., a double-headed or bidirectional arrow). Here, \code{A} and
        \code{B} are variable names in the model. If a name does not correspond
        to an observed variable, then it is assumed to be a latent variable.
        Spaces can appear freely in an arrow specification, and
        there can be any number of hyphens in the arrows, including zero: Thus,
        e.g., \code{A->B}, \code{A --> B}, and \code{A>B} are all legitimate
        and equivalent.}
      \item{2. Parameter name:}{The name of the regression coefficient, variance,
        or covariance specified by the arrow. Assigning the same name to two or
        more arrows results in an equality constraint. Specifying the parameter name
        as \code{NA} produces a fixed parameter.}
      \item{3. Value:}{start value for a free parameter or value of a fixed parameter.
        If given as \code{NA}, \code{sem} will compute the start value.}
      }
      
  Lines may end in a comment following \code{\#}.
      
  See \code{\link{sem}} for further details on model specification.
}

\value{
  An object of class \code{mod}, suitable as input for \code{\link{sem}}.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{
    \dontrun{
model.dhp <- specify.model()
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

model.dhp
    }
}

\keyword{models}
