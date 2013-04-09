\name{effects.sem}
\alias{effects.sem}
\alias{effects.msem}
\alias{print.semeffects}
\alias{print.semeffectsList}

\title{
Total, Direct, and Indirect Effects for Structural Equation Models
}

\description{
The \code{sem} method for the standard generic function \code{effects} computes total, direct, 
and indirect effects for a fitted structural equation model according to the method described in Fox (1980).
}

\usage{
\method{effects}{sem}(object, ...)
\method{effects}{msem}(object, ...)

\method{print}{semeffects}(x, digits = getOption("digits"), ...)
\method{print}{semeffectsList}(x, digits = getOption("digits"), ...)
}

\arguments{
  \item{object}{a fitted structural-equation model object produced by the \code{\link{sem}} function.}
  \item{x}{an object of class \code{semeffects} or \code{semeffectsList}, produced by \code{effects}.}
  \item{digits}{digits to print.}
  \item{\dots}{not used.}
}

\value{
\code{effect.sem} returns an object of class \code{semeffects} with \code{Total}, \code{Direct}, and \code{Indirect} elements.
}

\references{
  Fox, J. (1980)
  Effect analysis in structural equation models: Extensions and simplified methods of computation.
  \emph{Sociological Methods and Research}
  \bold{9}, 3--28.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{
    \dontrun{

# These examples are from Fox (1980)

# In the first pair of examples, readMoments() and specifyModel() read from the
# input stream. These examples cannot be executed via example() but can be entered
# at the command prompt. The Blau and Duncan example is repeated using file input;
# this example can be executed via example(). 

# The recursive Blau and Duncan basic stratification model:
#  x1 is father's education, x2 father's SES, y3 respondent's education,
#  y4 SES of respondent's first job, y5 respondent's SES in 1962

R.bd <- readMoments(names=c("x1", "x2", "y3", "y4", "y5"))
1
.516 1
.453 .438 1
.332 .417 .538 1
.322 .405 .596 .541 1

mod.bd <- specifyModel()
y3 <- x1, gam31
y3 <- x2, gam32
y4 <- x2, gam42
y4 <- y3, beta43
y5 <- x2, gam52
y5 <- y3, beta53
y5 <- y4, beta54

sem.bd <- sem(mod.bd, R.bd, N=20700, fixed.x=c("x1", "x2"))
summary(sem.bd)
effects(sem.bd)


# The nonrecursive Duncan, Haller, and Portes peer-influences model for observed variables:

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

model.dhp <- specifyModel()
RIQ      -> ROccAsp, gam51,  NA
RSES     -> ROccAsp, gam52,  NA
FSES     -> FOccAsp, gam63,  NA
FIQ      -> FOccAsp, gam64,  NA
FOccAsp  -> ROccAsp, beta56, NA
ROccAsp  -> FOccAsp, beta65, NA
ROccAsp <-> ROccAsp, ps55,   NA
FOccAsp <-> FOccAsp, ps66,   NA
ROccAsp <-> FOccAsp, ps56,   NA


# Note: The following generates a warning because not all of the variables
#       in the correlation matrix are used
sem.dhp <- sem(model.dhp, R.DHP, 329,
                fixed.x=c('RIQ', 'RSES', 'FSES', 'FIQ'))
summary(sem.dhp)
effects(sem.dhp)
    }
    
## the following example may be executed via example()

etc <- system.file(package="sem", "etc") # path to data and model files

# The recursive Blau and Duncan basic stratification model:
#  x1 is father's education, x2 father's SES, y3 respondent's education,
#  y4 SES of respondent's first job, y5 respondent's SES in 1962

(R.bd <- readMoments(file=file.path(etc, "R-Blau-Duncan.txt"),
					names=c("x1", "x2", "y3", "y4", "y5")))
(mod.bd <- specifyModel(file=file.path(etc, "model-Blau-Duncan.txt")))
sem.bd <- sem(mod.bd, R.bd, N=20700, fixed.x=c("x1", "x2"))
summary(sem.bd)
effects(sem.bd)
}


\keyword{models}

