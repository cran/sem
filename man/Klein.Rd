\name{Klein}
\alias{Klein}
\docType{data}
\title{Klein's Data on the U. S. Economy}
\usage{data(Klein)}
\description{
Data for Klein's (1950) simple econometric model of the U. S. economy.

The \code{Klein} data frame has 22 rows and 10 columns.

}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Year}{1921--1941}
    \item{C}{consumption.}
    \item{P}{private profits.}
    \item{Wp}{private wages.}
    \item{I}{investment.}
    \item{K.lag}{capital stock, lagged one year.}
    \item{X}{equilibrium demand.}
    \item{Wg}{government wages.}
    \item{G}{government non-wage spending.}
    \item{T}{indirect business taxes and net exports.}
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

Klein$P.lag <- c(NA, Klein$P[-22])
Klein$X.lag <- c(NA, Klein$X[-22])

summary(tsls(C ~ P + P.lag + I(Wp + Wg), 
    instruments=~1 + G + T + Wg + I(Year - 1931) + K.lag + P.lag + X.lag,
    data=Klein))
    
summary(tsls(I ~ P + P.lag + K.lag,
    instruments=~1 + G + T + Wg + I(Year - 1931) + K.lag + P.lag + X.lag,
    data=Klein))
    
summary(tsls(Wp ~ X + X.lag + I(Year - 1931),
    instruments=~1 + G + T + Wg + I(Year - 1931) + K.lag + P.lag + X.lag,
    data=Klein))

}
\keyword{datasets}
