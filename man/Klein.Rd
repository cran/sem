\name{Klein}
\alias{Klein}
\non_function{}
\title{Klein's Data on the U. S. Economy}
\usage{data(Klein)}
\description{
Data for Klein's (1950) simple econometric model of the U. S. economy.

The \code{Klein} data frame has 22 rows and 10 columns.

}
\format{
  This data frame contains the following columns:
  \describe{
    \item{year}{1921--1941}
    \item{c}{consumption.}
    \item{p}{private profits.}
    \item{wp}{private wages.}
    \item{i}{investment.}
    \item{k.lag}{capital stock, lagged one year.}
    \item{x}{equilibrium demand.}
    \item{wg}{government wages.}
    \item{g}{government non-wage spending.}
    \item{t}{indirect business taxes and net exports.}
  }
}

\source{
    Greene, W. H. (1993)
    \emph{Econometric Analysis, Second Edition.}
    Macmillan.
}
\references{
    Klein, L. (1950)
    \emph{Economic Fluctuations in the United States 1921--1941.}
    Wiley.
}
\examples{
data(Klein)
}
\keyword{datasets}
