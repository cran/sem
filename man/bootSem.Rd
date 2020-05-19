\name{bootSem}
\alias{bootSem}
\alias{bootSem.sem}
\alias{bootSem.msem}
\alias{print.bootsem}
\alias{summary.bootsem}
\alias{print.summary.bootsem}

\title{Bootstrap a Structural Equation Model}
\description{
Bootstraps a structural equation model in an \code{sem} object (as returned by
the \code{\link{sem}} function).
}

\usage{
bootSem(model, ...)

\method{bootSem}{sem}(model, R=100, Cov=cov, data=model$data, 
    max.failures=10, show.progress=TRUE, ...)

\method{bootSem}{msem}(model, R=100, Cov=cov, data=model$data, 
    max.failures=10, show.progress=TRUE, ...)

\method{print}{bootsem}(x, digits=getOption("digits"), ...)

\method{summary}{bootsem}(object,
    type=c("perc", "bca", "norm", "basic", "none"), level=0.95, ...)
}

\arguments{
  \item{model}{an \code{sem} or \code{msem} object, produced by the \code{sem} function.}
  \item{R}{the number of bootstrap replications; the default is 100, which should
    be enough for computing standard errors, but not confidence intervals
    (except for the normal-theory intervals).}
  \item{Cov}{a function to compute the input covariance or moment matrix; the default is
    \code{\link{cov}}. Use \code{\link[stats]{cor}} if the model is fit
    to the correlation matrix. The function \code{\link[polycor]{hetcor}} in the 
    \code{polycor} package will compute product-moment, polychoric, and
    polyserial correlations among mixed continuous and ordinal variables
    (see the first example below for an illustration).}
  \item{data}{in the case of a \code{sem} (i.e., single-group) \code{model}, a data set in a form suitable for \code{Cov};
    for example, for the default \code{Cov=cov}, \code{data} may be a numeric data frame or a numeric matrix.
    In the case of an \code{msem} (i.e., multi-group) \code{model}, a list of data sets (again in the appropriate form),
    one for each group; in this case, bootstrapping is done within each group, treating the groups as strata.
    Note that the original observations are required, not just the covariance 
    matrix of the observed variables in the model. The default is the data set stored in the
    \code{sem} object, which will be present only if the model was fit to a data set rather than to
    a covariance or moment matrix, and may not be in a form suitable for \code{Cov}.}
  \item{max.failures}{maximum number of consecutive convergence failures before \code{bootSem} gives up.}
  \item{show.progress}{display a text progress bar on the console tracing the bootstrap replications.}
  \item{x, object}{an object of class \code{bootsem}.}
  \item{digits}{controls the number of digits to print.}
  \item{type}{type of bootstrapped confidence intervals to compute; the
    default is \code{"perc"} (percentile); see 
    \code{\link[boot]{boot.ci}} for details.}
  \item{level}{level for confidence intervals; default is \code{0.95}.}
  \item{...}{in \code{bootSem}, arguments to be passed to 
    \code{\link{sem}}; otherwise ignored.} 
}

\details{
  \code{bootSem} implements the nonparametric bootstrap, assuming an
  independent random sample. Convergence failures in the bootstrap resamples
  are discarded (and a warning printed); more than \code{max.failures} consecutive convergence failures (default, 10)
  result in an error. You can use the \code{\link[boot]{boot}} function
  in the \pkg{boot} package for more complex sampling schemes and additional options.
  
  Bootstrapping is implemented by resampling the observations in
  \code{data}, recalculating the input covariance matrix with \code{Cov},
  and refitting the model with \code{\link{sem}}, using the
  parameter estimates from the original sample as start-values.
  
  \bold{Warning:} the bootstrapping process can be very time-consuming.
}

\value{
  \code{bootSem} returns an object of class \code{bootsem}, which inherits
    from class \code{boot}, supported by the \pkg{boot} package. The returned
    object contains the following components:
    
  \item{t0}{the estimated parameters in the model fit to the original data set.}
  \item{t}{a matrix containing the bootstrapped estimates, one bootstrap
    replication per row.}
  \item{data}{the data to which the model was fit.}
  \item{seed}{the value of \code{.Random.seed} when \code{bootSem} was called.}  
  \item{statistic}{the function used to produce the bootstrap replications;
    this is always the local function \code{refit} from \code{bootSem}.}
  \item{sim}{always set to \code{"ordinary"}; see the documentation for the
    \code{\link[boot]{boot}} function.}
  \item{stype}{always set to \code{"i"}; see the documentation for the
    \code{\link[boot]{boot}} function.}
  \item{call}{the call of the \code{bootSem} function.}
  \item{weights}{a vector of length equal to the number of observations \eqn{N},
    with entries \eqn{1/N}. For a multi-group model, the weights in group \eqn{j} are \eqn{1/N_j},
    the inverse of the number of observations in the group.}
  \item{strata}{a vector of length \eqn{N} containing the number of the stratum
    to which each observation belongs; for a single-group model, all entries
    are 1.} 
}

\references{
    Davison, A. C., and Hinkley, D. V. (1997)
    \emph{Bootstrap Methods and their Application.} Cambridge.
    
    Efron, B., and Tibshirani, R. J. (1993)
    \emph{An Introduction to the Bootstrap.} Chapman and Hall.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[boot]{boot}}, \code{\link{sem}}}

\examples{
    \donttest{

# A simple confirmatory factor-analysis model using polychoric correlations.
#  The polycor package is required for the hetcor function.

if (require(polycor)){

# The following function returns correlations computed by hetcor,
#   and is used for the bootstrapping.

hcor <- function(data) hetcor(data, std.err=FALSE)$correlations

model.cnes <- specifyModel(text="
F -> MBSA2, lam1
F -> MBSA7, lam2
F -> MBSA8, lam3
F -> MBSA9, lam4
F <-> F, NA, 1
MBSA2 <-> MBSA2, the1
MBSA7 <-> MBSA7, the2
MBSA8 <-> MBSA8, the3
MBSA9 <-> MBSA9, the4
")

R.cnes <- hcor(CNES)

sem.cnes <- sem(model.cnes, R.cnes, N=1529)
summary(sem.cnes)
}

#  Note: this can take a minute:

set.seed(12345) # for reproducibility

system.time(boot.cnes <- bootSem(sem.cnes, R=100, Cov=hcor, data=CNES))
summary(boot.cnes, type="norm")  
# cf., standard errors to those computed by summary(sem.cnes)
    }
    
    \dontrun{  # because of long execution time

# An example bootstrapping a multi-group model

mod.hs <- cfa(text="
spatial: visual, cubes, paper, flags
verbal: general, paragrap, sentence, wordc, wordm
memory: wordr, numberr, figurer, object, numberf, figurew
math: deduct, numeric, problemr, series, arithmet
")

mod.mg <- multigroupModel(mod.hs, groups=c("Female", "Male")) 

sem.mg <- sem(mod.mg, data=HS.data, group="Gender",
              formula = ~ visual + cubes + paper + flags +
              general + paragrap + sentence + wordc + wordm +
              wordr + numberr + figurer + object + numberf + figurew +
              deduct + numeric + problemr + series + arithmet
              )

# Note: this example can take several minutes or more;
#       you can decrease R if you just want to see how it works:

set.seed(12345) # for reproducibility

system.time(boot.mg <- bootSem(sem.mg, R=100))
summary(boot.mg, type="norm")
# cf., standard errors to those computed by summary(sem.mg)
    }
}

\keyword{htest}
\keyword{models}
