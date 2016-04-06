\name{specifyModel}
\alias{specifyModel}
\alias{specifyEquations}
\alias{cfa}
\alias{multigroupModel}
\alias{print.semmod}
\alias{print.semmodList}
\alias{classifyVariables}
\alias{removeRedundantPaths}
\alias{combineModels}
\alias{combineModels.semmod}
\alias{update.semmod}
\alias{edit.semmod}

\title{Specify a Structural Equation Model}
\description{
  Create the RAM specification of a structural equation model.
}
\usage{
specifyModel(file="", text, exog.variances=FALSE, endog.variances=TRUE, covs, 
	suffix="", quiet=FALSE)

specifyEquations(file="", text, ...)

cfa(file="", text, covs=paste(factors, collapse=","), 
    reference.indicators=TRUE, raw=FALSE, 
    subscript=c("name", "number"), ...)

multigroupModel(..., groups=names(models), allEqual=FALSE)

classifyVariables(model)

removeRedundantPaths(model, warn=TRUE)
\method{combineModels}{semmod}(..., warn=TRUE)
\method{update}{semmod}(object, file = "", text, ...)
\method{edit}{semmod}(name, ...)

\method{print}{semmod}(x, ...)
\method{print}{semmodList}(x, ...)
}

\arguments{
  \item{file}{The (quoted) file from which to read the model specification,
    including the path to the file if it is not in the current directory. If
    \code{""} (the default) and the \code{text} argument is not supplied,
    then the specification is read from the standard
    input stream, and is terminated by a blank line.}
  \item{text}{The model specification given as a character string, as an alternative
    to specifying the ]code{file} argument or reading the model specification from
    the input stream --- e.g., when the session is not interactive and there is no
    standard input.}
  \item{exog.variances}{If \code{TRUE} (the default is \code{FALSE}), free variance
  	parameters are added for the exogenous variables that lack them.}
  \item{endog.variances}{If \code{TRUE} (the default), free error-variance parameters
    are added for the endogenous variables that lack them.}
  \item{covs}{optional: a character vector of one or more elements, with each element
  	giving a string of variable names, separated by commas. Variances and covariances
  	among all variables in each such string are added to the model. For confirmatory
  	factor analysis models specified via \code{cfa}, \code{covs} defaults to all of
  	the factors in the model, thus specifying all variances and covariances among these factors.
  	\emph{Warning}: \code{covs="x1, x2"} and \code{covs=c("x1", "x2")} are \emph{not}
  	equivalent: \code{covs="x1, x2"} specifies the variance of \code{x1}, the variance
  	of \code{x2}, \emph{and} their covariance, while \code{covs=c("x1", "x2")} specifies
  	the variance of \code{x1} and the variance of \code{x2} \emph{but not} their covariance.}
  \item{suffix}{a character string (defaulting to an empty string) to be 
    appended to each parameter name; this can be convenient for specifying multiple-group models.}
  \item{reference.indicators}{if \code{FALSE}, the default, variances of factors are
    set to 1 by \code{cfa}; if \code{TRUE}, variances of factors are free parameters
    to estimate from the data, and instead the first factor loading for each factor is set
    to 1 to identify the model.}
  \item{raw}{if \code{TRUE} (the default is \code{FALSE}), a path from \code{Intercept} to each
    observed variable is added to the model, and the raw second moment for \code{Intercept} is
    fixed to \code{1}. The \code{sem} function should then be called with \code{raw=TRUE}, and either
    supplied with a data set (via the \code{data} argument) or a raw-moment matrix (via the \code{S}
    argument).}
  \item{subscript}{The \dQuote{subscripts} to be appended to \code{lam} to name factor-loading parameters,
    either \code{"name"} (the default) to use the names of observed variables, or \code{"number"} to
    number the parameters serially within each factor. Using \code{"number"} produces shorter parameter
    names.}
  \item{quiet}{if \code{FALSE}, the default, then the number of input lines is reported and
    a message is printed suggesting that \code{specifyEquations} or \code{cfa} be used.}
  \item{x, model, object, name}{An object of class \code{semmod} or \code{semmodList}, as produced by
    \code{specifyModel} or \code{multigroupModel}.}
  \item{warn}{print a warning if redundant paths are detected.}
  \item{...}{For \code{multigroupModel}, one or more optionally named arguments each of
     which is a \code{semmod} object produced, e.g., by \code{specifyModel}, \code{specifyEquations},
     or \code{cfa}; if only one such model is given, then it will be used for all groups defined by
     the \code{groups} argument. If parameters have the same name in different groups, then they will
     be constrained to be equal.
     For \code{specifyEquations} and \code{cfa}, arguments (such as \code{covs},
  	 in the case of \code{specifyEquations}) to be 
     passed to \code{specifyModel}; for \code{combineModels}, \code{sem} objects; 
     ignored in the \code{update} and \code{print} methods.}
  \item{groups}{a character vector of names for the groups in a multigroup model; taken by default from
  	 the names of the \code{...} arguments.}
  \item{allEqual}{if \code{FALSE} (the default), then if only one model object is given for
  	 a multigroup model, all corresponding parameters in the groups will be distinct; if \code{TRUE},
  	 all corresponding parameters will be constrained to be equal.}
}

\details{
  The principal functions for model specification are \code{specifyModel},
  to specify a model in RAM (path) format via single- and double-headed arrows;
  \code{specifyEquations}, to specify a model in equation format, which is then
  translated by the function into RAM format; and \code{cfa}, for compact
  specification of simple confirmatory factor analysis models.
  
  \code{specifyModel}:
  
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
      
  Lines may end in a comment following \code{#}.
  
  \code{specifyEquations}:

  For \code{specifyEquations}, each input line is either a regression equation or the specification
  of a variance or covariance. Regression equations are of the form
  
  \code{y = par1*x1 + par2*x2 + ... + park*xk}
  
  where \code{y} and the \code{x}s are variables in the model (either observed or latent),
  and the \code{par}s are parameters. If a parameter is given as a numeric value (e.g.,
  \code{1}) then it is treated as fixed. Note that no \dQuote{error} variable is included in
  the equation; \dQuote{error variances} are specified via either the \code{covs} argument,
  via \code{V(y) = par} (see immediately below), or are added automatically to the model
  when, as by default, \code{endog.variances=TRUE}. A regression equation may be split over more
  than one input by breaking at a \code{+}, so that \code{+} is either the last non-blank character
  on a line or the first non-blank character on the subsequent line.
  
  Variances are specified in the form \code{V(var) = par} and covariances in the form
  \code{C(var1, var2) = par}, where the \code{var}s are variables (observed or unobserved) in
  the model. The symbols \code{V} and \code{C} may be in either lower- or upper-case. If \code{par}
  is a numeric value (e.g., \code{1}) then it is treated as fixed. In conformity with the RAM model,
  a variance or covariance for an endogenous variable in the model is an \dQuote{error} variance or
  covariance.
  
  \emph{Warning}: If the \code{covs} argument to \code{specifyEquations} is used to specify
  variances and covariances, please be aware that
  \code{covs="x1, x2"} and \code{covs=c("x1", "x2")} are \emph{not}
  equivalent: \code{covs="x1, x2"} specifies the variance of \code{x1}, the variance
  of \code{x2}, \emph{and} their covariance, while \code{covs=c("x1", "x2")} specifies
  the variance of \code{x1} and the variance of \code{x2} \emph{but not} their covariance.
  
  To set a start value for a free parameter, enclose the numeric start value in parentheses after the
  parameter name, as \code{parameter(value)}.
  
  \code{cfa}:
  
  For \code{cfa}, each input line includes the names of the variables, separated by commas, 
  that load on the corresponding factor; the name of the factor is given optionally at the beginning
  of the line, followed by a colon. If necessary, the variables that load on a factor may be continued
  across two or more input lines; in this case, each such line but the last must end in a comma. A
  variable may load on more than one factor (as long as the resulting model is identified, of course),
  but each factor may appear in only one input line (or set of input lines, if the variable list
  is continued onto the next line). 
  
  Equality constraints for factor loadings can be set by using equal-signs (\code{=}) rather than commas 
  to separate observed variable names. For example, \code{fac1: x1=x2=x3, x4=x5} sets the loadings
  for \code{x1}, \code{x2}, and \code{x3} equal to each other, and the loadings for \code{x4} and \code{x5}
  equal to each other. 
  
  Equality constraints among error variances can similarly be specified by using \code{var:} or \code{variance:}
  at the beginning of a line (actually, any character string beginning with \code{var} will do, and thus
  no factor name may begin with the characters \code{var}). For example, \code{var: x1=x2=x3, x4=x5} sets the
  error variances for \code{x1}, \code{x2}, and \code{x3} equal to each other, and the 
  error variances for \code{x4} and \code{x5} equal to each other. There may be several lines beginning with
  \code{var:}.
  
  If the argument \code{reference.indicators=FALSE}, the default,
  \code{cfa} will fix the variance of each factor to 1, and by
  default include covariances (i.e., correlations) among all pairs of factors. Alternatively,
  if \code{reference.indicators=TRUE}, then the factor variances are free parameters to be estimated
  from the data, and the first loading for each factor is set to 1 to identify the model. These two
  approaches produce equivalent models, with the same fit to the data, but alternative parametrizations.
  Specifying the argument \code{covs=NULL} implicitly fixes the factor intercorrelations to 0.
  
  See \code{\link{sem}} and the examples for further details on model specification.
  
  Other Functions:
  
  \code{classifyVariables} classifies the variables in a model as endogenous or exogenous.
  
  \code{combineModels} and \code{removeRedundantPaths} take \code{semmod} objects as arguments and do what their names imply.
  
  The \code{file} input argument to the \code{update} method for \code{semmod} objects, which by default comes from
  standard input, is a set of update directives, one per line. There are five kinds of directives. In each case
  the directive begins with the directive name, followed by one or more fields separated by commas.
  
  \describe{
  
  	\item{1. delete:}{Remove a path from the model. Example: \code{delete, RSES -> FGenAsp}}
  	
  	\item{2. add:}{Add a path to the model. Example (the \code{NA} for the start value is optional): \code{add, RSES -> FGenAsp, gam14,  NA}}
  	
  	\item{3. replace:}{Replace every occurrence of the first string with the second in the variables and
  		parameters of the model. This directive may be used, for example, to change one variable to
  		another or to rename a parameter. Example: \code{replace, gam, gamma}, substitutes the string \code{"gamma"}
  		for \code{"gam"} wherever the latter appears, presumably in parameter names.}
  		
  	\item{4. fix:}{Fix a parameter that was formerly free. Example: \code{fix, RGenAsp  -> REdAsp, 1}}
  	
  	\item{5. free:}{Free a parameter that was formerly fixed. Example (the \code{NA} for the start value is optional): 
  	\code{free, RGenAsp  -> ROccAsp, lam11, NA}}
  }
  
  The \code{edit} method for \code{semmod} objects opens the model in the R editor.
}

\value{
  \code{specifyModel}, \code{specifyEquations}, \code{cfa}, \code{removeRedundantPaths}, \code{combineModels}, 
  \code{update}, and \code{edit} return an object of class \code{semmod}, suitable as input for \code{\link{sem}}.
  
  \code{multigroupModel} returns an object of class \code{semmodList}, also suitable as input for \code{\link{sem}}.
  
  \code{classifyVariables} returns a list with two character vectors: \code{endogenous}, containing the names of endogenous
  variables in the model; and \code{exogenous}, containing the names of exogenous variables.
}

\author{John Fox \email{jfox@mcmaster.ca} and Jarrett Byrnes}

\seealso{\code{\link{sem}}}

\examples{

# example using the text argument:

model.dhp <- specifyModel(text="
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
")    
model.dhp

   # same model in equation form:
model.dhp.1 <- specifyEquations(covs="RGenAsp, FGenAsp", text="
RGenAsp = gam11*RParAsp + gam12*RIQ + gam13*RSES + gam14*FSES + beta12*FGenAsp
FGenAsp = gam23*RSES + gam24*FSES + gam25*FIQ + gam26*FParAsp + beta21*RGenAsp
ROccAsp = 1*RGenAsp
REdAsp = lam21(1)*RGenAsp  # to illustrate setting start values
FOccAsp = 1*FGenAsp
FEdAsp = lam42(1)*FGenAsp
")
model.dhp

# Note: The following examples can't be run via example() because the 
#  default file argument requires that the model specification be entered
#  at the command prompt. The examples can be copied and run in an interactive 
#  session in the R console, however.

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

# A three-factor CFA model for the Thurstone mental-tests data, 
#    specified three equivalent ways:

R.thur <- readMoments(diag=FALSE, 
    names=c('Sentences','Vocabulary',
            'Sent.Completion','First.Letters','4.Letter.Words','Suffixes',
            'Letter.Series','Pedigrees', 'Letter.Group'))
.828                                              
.776   .779                                        
.439   .493    .46                                 
.432   .464    .425   .674                           
.447   .489    .443   .59    .541                    
.447   .432    .401   .381    .402   .288              
.541   .537    .534   .35    .367   .32   .555        
.38   .358    .359   .424    .446   .325   .598   .452

	#  (1a) in CFA format:

mod.cfa.thur.c <- cfa(reference.indicators=FALSE)
FA: Sentences, Vocabulary, Sent.Completion
FB: First.Letters, 4.Letter.Words, Suffixes
FC: Letter.Series, Pedigrees, Letter.Group

cfa.thur.c <- sem(mod.cfa.thur.c, R.thur, 213)
summary(cfa.thur.c)

	#  (1b) in CFA format, using reference indicators:
	
mod.cfa.thur.r <- cfa()
FA: Sentences, Vocabulary, Sent.Completion
FB: First.Letters, 4.Letter.Words, Suffixes
FC: Letter.Series, Pedigrees, Letter.Group

cfa.thur.r <- sem(mod.cfa.thur.r, R.thur, 213)
summary(cfa.thur.r)

	#  (2) in equation format:

mod.cfa.thur.e <- specifyEquations(covs="F1, F2, F3")
Sentences = lam11*F1
Vocabulary = lam21*F1
Sent.Completion = lam31*F1
First.Letters = lam42*F2
4.Letter.Words = lam52*F2
Suffixes = lam62*F2
Letter.Series = lam73*F3
Pedigrees = lam83*F3
Letter.Group = lam93*F3
V(F1) = 1
V(F2) = 1
V(F3) = 1

cfa.thur.e <- sem(mod.cfa.thur.e, R.thur, 213)
summary(cfa.thur.e)

	#  (3) in path format:

mod.cfa.thur.p <- specifyModel(covs="F1, F2, F3")
F1 -> Sentences,                      lam11
F1 -> Vocabulary,                     lam21
F1 -> Sent.Completion,                lam31
F2 -> First.Letters,                  lam41
F2 -> 4.Letter.Words,                 lam52
F2 -> Suffixes,                       lam62
F3 -> Letter.Series,                  lam73
F3 -> Pedigrees,                      lam83
F3 -> Letter.Group,                   lam93
F1 <-> F1,                            NA,     1
F2 <-> F2,                            NA,     1
F3 <-> F3,                            NA,     1

cfa.thur.p <- sem(mod.cfa.thur.p, R.thur, 213)
summary(cfa.thur.p)

# The Thursstone CFA model with equality constraints on the
#  factor loadings and error variances

mod.cfa.thur.ceq <- cfa(reference.indicators=FALSE)
FA: Sentences = Vocabulary = Sent.Completion
FB: First.Letters = 4.Letter.Words = Suffixes
FC: Letter.Series = Pedigrees = Letter.Group
var: Sentences = Vocabulary = Sent.Completion
var: First.Letters = 4.Letter.Words = Suffixes
var: Letter.Series = Pedigrees = Letter.Group

cfa.thur.ceq <- sem(mod.cfa.thur.ceq, R.thur, 213)
summary(cfa.thur.ceq)
anova(cfa.thur.c, cfa.thur.ceq)
pathDiagram(cfa.thur.ceq, ignore.double=FALSE, ignore.self=TRUE,
    min.rank="FA, FB, FC", edge.labels="values")

# a multigroup CFA model fit to the Holzinger-Swineford
#   mental-tests data (from the MBESS package)

library(MBESS)
data(HS.data)

mod.hs <- cfa()
spatial: visual, cubes, paper, flags
verbal: general, paragrap, sentence, wordc, wordm
memory: wordr, numberr, figurer, object, numberf, figurew
math: deduct, numeric, problemr, series, arithmet

mod.mg <- multigroupModel(mod.hs, groups=c("Female", "Male")) 

sem.mg <- sem(mod.mg, data=HS.data, group="Gender",
              formula = ~ visual + cubes + paper + flags +
              general + paragrap + sentence + wordc + wordm +
              wordr + numberr + figurer + object + numberf + figurew +
              deduct + numeric + problemr + series + arithmet
              )
summary(sem.mg)

	# with cross-group equality constraints:
	
mod.mg.eq <- multigroupModel(mod.hs, groups=c("Female", "Male"), allEqual=TRUE)

sem.mg.eq <- sem(mod.mg.eq, data=HS.data, group="Gender",
              formula = ~ visual + cubes + paper + flags +
                general + paragrap + sentence + wordc + wordm +
                wordr + numberr + figurer + object + numberf + figurew +
                deduct + numeric + problemr + series + arithmet
              )
summary(sem.mg.eq)
    }
}

\keyword{models}
