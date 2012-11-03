\name{miSem}
\alias{miSem}
\alias{miSem.semmod}
\alias{miSem.semmodList}
\alias{print.miSem}
\alias{summary.miSem}

\title{
Estimate a Structural Equation Model By Multiple Imputation
}

\description{
\code{miSem} uses the \code{\link[mi]{mi}} function in the \pkg{mi} package to generate multiple imputations of missing
data, fitting the specified model to each completed data set.
}

\usage{
miSem(model, ...)

\method{miSem}{semmod}(model, ..., data, formula = ~., raw = FALSE, 
    fixed.x=NULL, objective=objectiveML, n.imp=5, n.iter=30, 
    seed=sample(1e+06, 1), mi.args=list())
    
\method{miSem}{semmodList}(model, ..., data, formula = ~., group, raw=FALSE, 
        fixed.x=NULL, objective=msemObjectiveML,
        n.imp=5, n.iter=30, seed=sample(1e6, 1), mi.args=list())

\method{print}{miSem}(x, ...)

\method{summary}{miSem}(object, digits=max(3, getOption("digits") - 2), ...)
}

\arguments{
  \item{model}{
    an SEM model-description object of class \code{semmod} or \code{semmodList}, created by \code{\link{specifyEquations}}
    \code{\link{cfa}}, or \code{\link{specifyModel}}, in the case of a multi-group model 
    in combination with \code{\link{multigroupModel}}.
    }
  \item{\dots, formula, raw, fixed.x, objective, group}{
    arguments to be passed to \code{\link{sem}}.
    }
  \item{data}{
    an R data frame, presumably with some missing data (encoded as \code{NA}), containing the data for fitting
    the SEM, possibly along with other variables to use to obtain multiple imputations of missing values. In
    the case of a multi-group model, this must be a \emph{single} data frame.
    }

  \item{n.imp}{
    number of imputations (default \code{5}).
    }

  \item{n.iter}{
    number of iterations for the multiple-imputation process (default \code{30}).
    }

  \item{seed}{
    seed for the random-number generator (default is an integer sampled from 1 to 1E6); stored in the resulting object.
    }

  \item{mi.args}{
    other arguments to be passed to \code{\link[mi]{mi}}.
    }
  \item{x, object}{
    an object of class \code{"miSem"}.
    }
  \item{digits}{
    for printing numbers.
  }
}

\value{
    \code{miSem} returns an object of class \code{"miSem"} with the following components:
    \item{initial.fit}{an \code{sem} model object produced using \code{\link{objectiveFIML}} if \code{raw=TRUE},
    or the objective function given by the \code{objective} argument otherwise.}
    \item{mi.fits}{a list of \code{sem} model objects, one for each imputed data set.}
    \item{imputation}{the object produced by \code{\link[mi]{mi}}, containing the imputed data sets.}
    \item{seed}{the seed used for the random number generator.}
}

\references{
    Yu-Sung Su, Andrew Gelman, Jennifer Hill, Masanao Yajima. (2011). 
    ``Multiple imputation with diagnostics (mi) in R: Opening windows into the black box.'' 
    \emph{Journal of Statistical Software} 45(2).
}
\author{
    John Fox \email{jfox@mcmaster.ca}
    }

\seealso{
\code{\link{sem}}, \code{\link[mi]{mi}}
}
\examples{
    \dontrun{
mod.cfa.tests <- cfa(raw=TRUE)
verbal: x1, x2, x3
math: y1, y2, y3

imps <- miSem(mod.cfa.tests, data=Tests, fixed.x="Intercept", raw=TRUE, seed=12345,
              mi.args=list(add.noise=noise.control(post.run.iter=30)))
summary(imps, digits=3) 

library(MBESS) # for data
data(HS.data)

# introduce some missing data:
HS <- HS.data[, c(2,7:10,11:15,20:25,26:30)]
set.seed(12345)
r <- sample(301, 100, replace=TRUE)
c <- sample(2:21, 100, replace=TRUE)
for (i in 1:100) HS[r[i], c[i]] <- NA

mod.hs <- cfa()
spatial: visual, cubes, paper, flags
verbal: general, paragrap, sentence, wordc, wordm
memory: wordr, numberr, figurer, object, numberf, figurew
math: deduct, numeric, problemr, series, arithmet

mod.mg <- multigroupModel(mod.hs, groups=c("Female", "Male")) 
imps.mg <- miSem(mod.mg, data=HS, group="Gender",
             seed=12345, n.iter=50, 
             mi.args=list(add.noise=noise.control(post.run.iter=30)))
summary(imps.mg, digits=3)
    }
}

\keyword{models}
