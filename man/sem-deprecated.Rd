\name{sem-deprecated}
\alias{sem-deprecated}

\alias{boot.sem}
\alias{mod.indices}
\alias{normalized.residuals}
\alias{path.diagram}
\alias{raw.moments}
\alias{read.moments}
\alias{specify.model}
\alias{standardized.coefficients}
\alias{standardized.residuals}
\alias{std.coef}

\title{Deprecated Functions in the sem Package}

\description{
  These functions are provided for compatibility with older versions of the \pkg{sem} package only, 
  and may be removed eventually. Although an effort has been made to insure backwards-compatibility,
  commands that worked in versions of the \pkg{sem} package prior to version 2.0-0
  will not necessarily work in version 2.0-0 and beyond, or may not work in the same manner.
}

\usage{
boot.sem(...)
mod.indices(...)
normalized.residuals(...)
path.diagram(...)
raw.moments(...)
read.moments(...)
specify.model(...)
standardized.coefficients(...)
standardized.residuals(...)
std.coef(...)
}

\arguments{
  \item{\dots}{pass arguments down to replacements for deprecated functions.}
}

\details{

\code{boot.sem} is now a synonym for the \code{\link{bootSem}} function.

\code{mod.indices} is now a synonym for \code{\link{modIndices}}.

\code{normalized.residuals} is now a synonym for \code{\link{normalizedResiduals}}.

\code{path.diagram} is now a synonym for \code{\link{pathDiagram}}.

\code{raw.moments} is now a synonym for \code{\link{rawMoments}}.

\code{read.moments} is now a synonym for \code{\link{readMoments}}.

\code{specify.model} is now a synonym for \code{\link{specifyModel}}.

\code{standardized.coefficients} and \code{std.coef} are now synonyms for the \code{\link{standardizedCoefficients}} and \code{\link{stdCoef}} 
functions.

\code{standardized.residuals} is now a synonym for \code{\link{standardizedResiduals}}.

}