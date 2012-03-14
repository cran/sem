\name{ram}
\alias{ram}

\title{RAM Matrix for a Structural-Equation Model}
\description{
  Print the labelled RAM definition matrix for a structural-equation
  model fit by \code{sem}.
}

\usage{
ram(object, digits=5, startvalues=FALSE)
}


\arguments{
  \item{object}{an object of class \code{sem} returned by the \code{sem} function.}
  \item{digits}{number of digits for printed output.}
  \item{startvalues}{if \code{TRUE}, start values for parameters are printed;
    otherwise, the parameter estimates are printed; the default is \code{FALSE}.}
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
    }
}
\keyword{models}
