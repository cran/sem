\name{information.criteria}
\alias{information.criteria}
\alias{AICc}
\alias{CAIC}

\title{
Additional Information Criteria
}

\description{
These are generic functions for computing, respectively, the AICc (second-order corrected Akaike Information Criterion) and CAIC (consistent Akaike Information Criterion).
}

\usage{
AICc(object, ...)

CAIC(object, ...)
}

\arguments{
  \item{object}{an object for which an appropriate \code{AICc} or \code{CAIC} method exists.}
  \item{\ldots}{possible additional arguments for methods.}
}

\author{Jarrett Byrnes and John Fox \email{jfox@mcmaster.ca}}

\references{
     Burnham, K. P., and Anderson, D. R. (1998)
    \emph{Model Selection and Inference: A Practical Information-Theoretical Approach.} Springer.
    
     Bozdogan, H. (1987)
     Model selection and Akaike's information criterion (AIC).
	\emph{Psychometrika}  bold{52}, 345--370.
}

\seealso{
\code{\link{AICc.objectiveML}}, \code{\link{CAIC.objectiveML}}
}

\keyword{models}