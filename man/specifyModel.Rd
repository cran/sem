\name{specifyModel}
\alias{specifyModel}
\alias{specifyEquations}
\alias{print.semmod}
\alias{classifyVariables}
\alias{removeRedundantPaths}
\alias{combineModels}
\alias{combineModels.semmod}
\alias{update.semmod}

\title{Specify a Structural Equation Model}
\description{
  Create the RAM specification of a structural equation model.
}
\usage{
specifyModel(file="", exog.variances=FALSE, endog.variances=TRUE, covs, quiet=FALSE)

specifyEquations(file="", ...)

classifyVariables(model)

removeRedundantPaths(model, warn=TRUE)
\method{combineModels}{semmod}(..., warn=TRUE)
\method{update}{semmod}(object, file = "", ...)

\method{print}{semmod}(x, ...)
}

\arguments{
  \item{file}{The (quoted) file from which to read the model specification,
    including the path to the file if it is not in the current directory. If
    \code{""} (the default), then the specification is read from the standard
    input stream, and is terminated by a blank line.}
  \item{exog.variances}{If \code{TRUE} (the default is \code{FALSE}), free variance
  	parameters are added for the exogenous variables that lack them.}
  \item{endog.variances}{If \code{TRUE} (the default), free error-variance parameters
    are added for the endogenous variables that lack them.}
  \item{covs}{optional: a character vector of one or more elements, with each element
  	giving a string of variable names, separated by commas. Variances and covariances
  	among all variables in each such string are added to the model.}
  \item{quiet}{if \code{FALSE}, the default, then the number of input lines is reported.}
  \item{x, model, object}{An object of class \code{semmod}, as produced by
    \code{specifyModel}.}
  \item{warn}{print a warning if redundant paths are detected.}
  \item{...}{for \code{specifyEquations}, arguments such as \code{covs} to be 
     passed to \code{specifyModel}; for \code{combineModels}, \code{sem} objects; 
     ignored in the \code{update} and \code{print} methods.}
}

\details{
  The principal functions for model specification are \code{specifyModel},
  to specify a model in RAM (path) format via single- and double-headed arrows; and
  \code{specifyEquations}, to specify a model in equation format, which is then
  translated by the function into RAM format.
  
  Each line of the RAM specification for \code{specifyModel} consists of three (unquoted) entries, 
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
        If given as \code{NA} (or simply omitted), \code{sem} will compute the start value.}
      }
      
  Lines may end in a comment following \code{\#}.

  For \code{specifyEquations}, each input line is either a regression equation or the specification
  of a variance or covariance. Regression equations are of the form
  
  \code{y = par1*x1 + par2*x2 + ... + park*xk}
  
  where \code{y} and the \code{x}s are variables in the model (either observed or latent),
  and the \code{par}s are parameters. If a parameter is given as a numeric value (e.g.,
  \code{1}) then it is treated as fixed. Note that no \dQuote{error} variable is included in
  the equation; \dQuote{error variances} are specified via either the \code{covs} argument,
  via \code{V(y) = par} (see immediately below), or are added automatically to the model
  when, as by default, \code{endog.variances=TRUE}.
  
  Variances are specified in the form \code{V(var) = par} and covariances in the form
  \code{C(var1, var2) = par}, where the \code{var}s are variables (observed or unobserved) in
  the model. The symbols \code{V} and \code{C} may be in either lower- or upper-case. If \code{par}
  is a numeric value (e.g., \code{1}) then it is treated as fixed. In conformity with the RAM model,
  a variance or covariance for an endogenous variable in the model is an \dQuote{error} variance or
  covariance.
  
  To set a start value for a free parameter, enclose the numeric start value in parentheses after the
  parameter name, as \code{parameter(value)}.
  
  See \code{\link{sem}} for further details on model specification.
  
  \code{classifyVariables} classifies the variables in a model as endogenous or exogenous.
  
  \code{combineModels} and \code{removeRedundantPaths} take \code{semmod} objects as arguments and do what their names imply.
  
  The \code{file} input argument to the \code{update} method for \code{semmod} objects, which by default comes from
  standard input, is a set of update directives, one per line. There are three kinds of directives. In each case
  the directive begins with the directive name, followed by one or more fields separated by commas.
  
  \describe{
  
  	\item{1. delete:}{Remove a path from the model. Example: \code{delete, RSES -> FGenAsp}}
  	
  	\item{2. add:}{Add a path to the model. Example: \code{add, RSES -> FGenAsp, gam14,  NA}}
  	
  	\item{3. replace:}{Replace every occurrence of the first string with the second in the variables and
  		parameters of the model. This directive may be used, for example, to change one variable to
  		another or to rename a parameter. Example: \code{replace, gam, gamma}, substitutes the string \code{"gamma"}
  		for \code{"gam"} wherever the latter appears, presumably in parameter names.}
  }
  
}

\value{
  \code{specifyModel}, \code{removeRedundantPaths}, \code{combineModels}, and \code{update} return an object of class \code{semmod}, suitable as input for \code{\link{sem}}.
  
  \code{classifyVariables} returns a list with two character vectors: \code{endogenous}, containing the names of endogenous
  variables in the model; and \code{exogenous}, containing the names of exogenous variables.
}

\author{John Fox \email{jfox@mcmaster.ca} and Jarrett Byrnes}

\seealso{\code{\link{sem}}}

\examples{

# Note: These examples can't be run via example() because the default file
#  argument of specifyModel() requires that the model specification be entered
#  at the command prompt. The examples can be copied and run in the R console,
#  however.

    \dontrun{
model.dhp <- specifyModel()
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
    
# an equivalent specification, allowing specifyModel() to generate
#  variance parameters for endogenous variables (and suppressing
#  the unnecessary trailing NAs):
 
model.dhp <- specifyModel()
RParAsp  -> RGenAsp, gam11
RIQ      -> RGenAsp, gam12
RSES     -> RGenAsp, gam13
FSES     -> RGenAsp, gam14
RSES     -> FGenAsp, gam23
FSES     -> FGenAsp, gam24
FIQ      -> FGenAsp, gam25
FParAsp  -> FGenAsp, gam26
FGenAsp  -> RGenAsp, beta12
RGenAsp  -> FGenAsp, beta21
RGenAsp  -> ROccAsp,  NA,     1
RGenAsp  -> REdAsp,  lam21
FGenAsp  -> FOccAsp,  NA,     1
FGenAsp  -> FEdAsp,  lam42
RGenAsp <-> FGenAsp, ps12

model.dhp

# Another equivalent specification, telling specifyModel to add paths for 
#   variances and covariance of RGenAsp and FGenAsp:
 
model.dhp <- specifyModel(covs="RGenAsp, FGenAsp")
RParAsp  -> RGenAsp, gam11
RIQ      -> RGenAsp, gam12
RSES     -> RGenAsp, gam13
FSES     -> RGenAsp, gam14
RSES     -> FGenAsp, gam23
FSES     -> FGenAsp, gam24
FIQ      -> FGenAsp, gam25
FParAsp  -> FGenAsp, gam26
FGenAsp  -> RGenAsp, beta12
RGenAsp  -> FGenAsp, beta21
RGenAsp  -> ROccAsp,  NA,     1
RGenAsp  -> REdAsp,  lam21
FGenAsp  -> FOccAsp,  NA,     1
FGenAsp  -> FEdAsp,  lam42

model.dhp

# The same model in equation format:

model.dhp.1 <- specifyEquations(covs="RGenAsp, FGenAsp")
RGenAsp = gam11*RParAsp + gam12*RIQ + gam13*RSES + gam14*FSES + beta12*FGenAsp
FGenAsp = gam23*RSES + gam24*FSES + gam25*FIQ + gam26*FParAsp + beta21*RGenAsp
ROccAsp = 1*RGenAsp
REdAsp = lam21(1)*RGenAsp  # to illustrate setting start values
FOccAsp = 1*FGenAsp
FEdAsp = lam42(1)*FGenAsp

model.dhp

classifyVariables(model.dhp)

# updating the model to impose equality constraints
#  and to rename the latent variables and gamma parameters

model.dhp.eq <- update(model.dhp)
delete, RSES -> FGenAsp
delete, FSES -> FGenAsp
delete, FIQ  -> FGenAsp
delete, FParAsp -> FGenAs
delete, RGenAsp  -> FGenAsp
add, RSES     -> FGenAsp, gam14,  NA
add, FSES     -> FGenAsp, gam13,  NA
add, FIQ      -> FGenAsp, gam12,  NA
add, FParAsp  -> FGenAsp, gam26,  NA
add, RGenAsp  -> FGenAsp, beta12, NA
replace, gam, gamma
replace, Gen, General

model.dhp.eq

    }
}

\keyword{models}
