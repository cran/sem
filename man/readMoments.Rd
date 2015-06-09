\name{readMoments}
\alias{readMoments}

\title{Input a Covariance, Correlation, or Raw Moment Matrix}
    
\description{
This functions makes it simpler to input covariance, correlation, and raw-moment 
matrices to be analyzed by the \code{\link{sem}} function. The matrix
is input in lower-triangular form on as many lines as is convenient, 
omitting the above-diagonal elements. The
elements on the diagonal may also optionally be omitted, in which case they
are taken to be 1.
}

\usage{
readMoments(file="", text, diag=TRUE, 
    names=as.character(paste("X", 1:n, sep = "")))
}

\arguments{
  \item{file}{The (quoted) file from which to read the moment matrix,
    including the path to the file if it is not in the current directory. If
    \code{""} (the default) and the \code{text} argument is absent, 
    then the moment matrix is read from the standard
    input stream, and is terminated by a blank line.}
  \item{text}{The moment matrix given as a character string, as an alternative
    to specifying the \code{file} argument or reading the moments from
    the input stream --- e.g., when the session is not interactive and there is no
    standard input.}
  \item{diag}{If \code{TRUE} (the default), then the input matrix includes
    diagonal elements.}
  \item{names}{a character vector containing the names of the variables, to
    label the rows and columns of the moment matrix.}
}

\value{
 Returns a lower-triangular matrix (i.e., with zeroes above the main diagonal)
 suitable for input to \code{sem}.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{

R.DHP <- readMoments(diag=FALSE, names=c("ROccAsp", "REdAsp", "FOccAsp", 
                "FEdAsp", "RParAsp", "RIQ", "RSES", "FSES", "FIQ", "FParAsp"),
                text="
    .6247     
    .3269  .3669       
    .4216  .3275  .6404
    .2137  .2742  .1124  .0839
    .4105  .4043  .2903  .2598  .1839
    .3240  .4047  .3054  .2786  .0489  .2220
    .2930  .2407  .4105  .3607  .0186  .1861  .2707
    .2995  .2863  .5191  .5007  .0782  .3355  .2302  .2950
    .0760  .0702  .2784  .1988  .1147  .1021  .0931 -.0438  .2087
")   
R.DHP

#the following will work only in an interactive sessions:
    \dontrun{
R.DHP <- readMoments(diag=FALSE, names=c("ROccAsp", "REdAsp", "FOccAsp", 
                "FEdAsp", "RParAsp", "RIQ", "RSES", "FSES", "FIQ", "FParAsp"))
    .6247     
    .3269  .3669       
    .4216  .3275  .6404
    .2137  .2742  .1124  .0839
    .4105  .4043  .2903  .2598  .1839
    .3240  .4047  .3054  .2786  .0489  .2220
    .2930  .2407  .4105  .3607  .0186  .1861  .2707
    .2995  .2863  .5191  .5007  .0782  .3355  .2302  .2950
    .0760  .0702  .2784  .1988  .1147  .1021  .0931 -.0438  .2087
    
R.DHP
    }
}

\keyword{manip}

