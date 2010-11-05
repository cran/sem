\name{sem} 

\alias{sem}
\alias{sem.semmod}
\alias{sem.default}
\alias{startvalues}
\alias{coef.sem}
\alias{vcov.sem}
\alias{df.residual.sem}

\title{General Structural Equation Models}
\description{
  \code{sem} fits general structural equation models (with both observed and
  unobserved variables). Observed variables are also called \emph{indicators} or 
  \emph{manifest variables}; unobserved variables are also called \emph{factors}
  or \emph{latent variables}. Normally, the generic function (\code{sem}) is
  called directly with a \code{semmod} first argument, invoking the \code{sem.semmod}
  method, which in turn sets up a call to the \code{sem.default} method; thus, the user
  may wish to specify arguments accepted by the \code{semmod} and \code{default} methods.
}

\usage{
\method{sem}{semmod}(model, S, N, data, raw=FALSE, obs.variables=rownames(S), 
		fixed.x=NULL, formula= ~ ., robust=!missing(data), debug=FALSE, ...)
    
\method{sem}{default}(model, S, N, data=NULL, raw=FALSE, param.names, 
		var.names, fixed.x=NULL, robust=!is.null(data), semmod=NULL, debug=FALSE,
		analytic.gradient=TRUE, warn=FALSE, maxiter=500, 
		par.size=c("ones", "startvalues"), start.tol=1E-6, 
		optimizer=optimizerNlm, objective=objectiveML, ...) 
    
startvalues(S, ram, debug=FALSE, tol=1E-6)

\method{coef}{sem}(object, standardized=FALSE, ...)
\method{vcov}{sem}(object, robust=FALSE, analytic=inherits(object, "objectiveML"), ...)
\method{df.residual}{sem}(object, ...)
}

\arguments{
  \item{model}{RAM specification, which is a simple encoding of the path
    diagram for the model. The model may be given either in symbolic 
    form (as a \code{semmod} object, as returned by the \code{\link{specifyModel}} function,
    or as a character matrix), invoking \code{sem.semmod}, which calls \code{sem.default} 
    after setting up the model,
    or (less conveniently) in numeric form, invoking \code{sem.default} directly, which is
    not recommended (see \bold{Details} below).}
  \item{S}{covariance matrix among observed variables; may be input as a symmetric matrix,
    or as a lower- or upper-triangular matrix. \code{S} may also be a raw (i.e., \dQuote{uncorrected})
    moment matrix --- that is, a sum-of-squares-and-products matrix divided by \code{N}. This
    form of input is useful for fitting models with intercepts, in which case the moment matrix
    should include the mean square and cross-products for a unit variable all of whose entries are 1;
    of course, the raw mean square for the unit variable is 1. Raw-moment matrices may be computed
    by \code{\link{rawMoments}}. If the \code{ram} argument is given in symbolic form, then
    the observed-variable covariance or raw-moment matrix may contain variables that do not appear in the model,
    in which case a warning is printed. As an alternative to specifying \code{S} the user may
    supply a data frame containing the data for the model (see the argument \code{data}).}
  \item{N}{number of observations on which the covariance matrix is based.}
  \item{data}{As a generally preferable alternative to specifying \code{S} and 
 	\code{N}, the user may supply a data frame containing the data to which the model is to be fit.}
  \item{raw}{\code{TRUE} if \code{S} is a raw moment matrix or if a raw moment matrix ---  as
    opposed to a covariance matrix --- is to be
  	computed from \code{data}; the default is \code{FALSE}.}
  \item{obs.variables}{names of observed variables, by default taken from the row names of
    the covariance or moment matrix \code{S}, which may be given directly or generated according to the
    \code{data} and \code{formula} arguments.}
  \item{fixed.x}{names (if the \code{ram} matrix is given in symbolic form) or indices 
    (if it is in numeric form) of fixed exogenous variables. Specifying these obviates
    the necessity of having to fix the variances and covariances among these 
    variables (and produces correct degrees of freedom for the model chisquare).}
  \item{formula}{a one-sided formula, to be applied to \code{data} to generate the variables for
  	which covariances or raw moments are computed. The default formula is \code{~.}, i.e., all of the
  	variables in the data, including an implied intercept; if a covariance matrix is to be computed,
  	the constant is suppressed.}
  \item{robust}{In \code{sem}: if \code{TRUE}, then quantities are calculated that can be used to compute robust
  	estimates of coefficient standard errors and robust tests when the model is fit by multinormal maximum likelihood;
  	the default is \code{TRUE} when the \code{data} argument is \code{TRUE}, and this option is only
  	available when the \code{data} argument is given. In \code{vcov}: if \code{TRUE}, return a robust coefficient covariance matrix
  	(if \code{object} contains the requisite information).}
  \item{semmod}{a \code{semmod} object containing the description of the model; optional, and normally supplied
  	not directly by the user but via the \code{semmod} method for \code{sem}.}
  \item{debug}{if \code{TRUE}, some information is printed to help you debug the symbolic
    model specification; for example, if a variable name is misspelled, \code{sem} will
    assume that the variable is a (new) latent variable. Information about
    the optimization will also be printed, but details will vary with the optimizer employed. 
     The default is \code{FALSE}.}
  \item{...}{arguments to be passed down, including from \code{sem.default} to the 
    optimizer.}
  \item{param.names}{names of the \eqn{t} free parameters, given in their numerical order;
    default names are \code{Param1}, \ldots, \code{Param}\emph{t}. Note: Should not be
    specified when the model is given in symbolic form.}
  \item{var.names}{names of the \eqn{m} entries of the \eqn{v} vector
    (typically the observed and latent variables --- see below), given in their
    numerical order; default names are \code{Var1}, \ldots, \code{Var}\emph{m}.
    Note: Should not be specified when the model is given in symbolic form.}
  \item{analytic.gradient}{if \code{TRUE} (the default), then analytic first derivatives are
    used in the maximization of the likelihood if the optimzer employed will accept them; 
    otherwise numeric derivatives are used, again if the optimizer will compute them.}
  \item{warn}{if \code{TRUE}, warnings produced by the optimization function will be printed.
    This should generally not be necessary, since \code{sem} prints its own warning, and saves
    information about convergence. The default is \code{FALSE}.}
  \item{maxiter}{the maximum number of iterations for the optimization of the objective function,
    to be passed to the optimizer.}
  \item{par.size}{the anticipated size of the free parameters; if \code{"ones"},
    a vector of ones is used; if \code{"startvalues"}, taken from the start values.
    You can try changing this argument if you encounter convergence problems.
    The default is \code{"startvalues"} if the largest input variance is at
    least 100 times the smallest, and \code{"ones"} otherwise. Whether this argument is actually
    used depends upon the optimizer employed.}
  \item{start.tol, tol}{if the magnitude of an automatic start value is less than \code{start.tol}, then
    it is set to \code{start.tol}; defaults to 1E-6.}
  \item{optimizer}{a function to be used to minimize the objective function; the default is 
  	\code{\link{optimizerNlm}}, which employs the standard R \code{\link[stats]{nlm}} optimizer.
  	Alternatives are \code{\link{optimizerOptim}}, which employs \code{\link[stats]{optim}},
  	and \code{\link{optimizerNlminb}}, which uses \code{\link[stats]{nlminb}} --- or
  	the user can supply an optimizer.}
  \item{objective}{An objective function to be minimized, sometimes called a \dQuote{fit} function 
  	in the SEM literature. The default is \code{\link{objectiveML}},
  	which produces maximum-likelihood estimates assuming multinormality. An alternative is 
  	\code{\link{objectiveGLS}}, which produced generalized least squares estimates --- or the user
  	can supply an objective function to be minimzed.}
  \item{ram}{numeric RAM matrix.}
  \item{object}{an object of class \code{"sem"}, returned by \code{sem}.}
  \item{standardized}{if \code{TRUE}, return standardized coefficients.}
  \item{analytic}{return an analytic (as opposed to numeric) estimate of the coefficient covariance matrix;
  	at present only available for the \code{\link{objectiveML}} objective function.}
} 

\details{
  The model is set up using RAM (\dQuote{reticular action model} -- don't ask!) 
  notation -- a simple format
  for specifying general structural equation models by coding the
  \dQuote{arrows} in the path diagram for the model (see, e.g., McArdle and McDonald, 1984).
  
  The variables in the \eqn{v} vector in the model (typically, the observed and
  unobserved variables, but not error variables) are numbered from 1 to \eqn{m}.  
  the RAM matrix contains one row for each (free or constrained) parameter of the model, and
  may be specified either in symbolic format or in numeric format.
  
  A symbolic \code{ram} matrix consists of three columns, as follows:
  
    \describe{
      \item{1. Arrow specification:}{This is a simple formula, of the form
        \code{"A -> B"} or, equivalently, \code{"B <- A"} for a regression
        coefficient (i.e., a single-headed or directional arrow); 
        \code{"A <-> A"} for a variance or \code{"A <-> B"} for a covariance
        (i.e., a double-headed or bidirectional arrow). Here, \code{A} and
        \code{B} are variable names in the model. If a name does not correspond
        to an observed variable, then it is assumed to be a latent variable.
        Spaces can appear freely in an arrow specification, and
        there can be any number of hyphens in the arrows, including zero: Thus,
        e.g., \code{"A->B"}, \code{"A --> B"}, and \code{"A>B"} are all legitimate
        and equivalent.}
      \item{2. Parameter name:}{The name of the regression coefficient, variance,
        or covariance specified by the arrow. Assigning the same name to two or
        more arrows results in an equality constraint. Specifying the parameter name
        as \code{NA} produces a fixed parameter.}
      \item{3. Value:}{start value for a free parameter or value of a fixed parameter.
        If given as \code{NA}, \code{sem} will compute the start value.}
      }
      
  It is simplest to construct the RAM matrix with the \code{\link{specifyModel}} function,
  which returns an object of class \code{semmod}, and also incorporates some model-specification
  convenience shortcuts. This process is illustrated in the examples below. 
  
  A numeric \code{ram} matrix consists of five columns, as follows:
    
    \describe{
        \item{1. Number of arrow heads:}{1 (directed arrow) or 2 (covariance).}
        \item{2. Arrow \emph{to}:}{index of the variable at the head of
            a directional arrow, or at one end of a bidirectional arrow.
            Observed variables should be assigned the numbers 1 to \eqn{n}, where
            \eqn{n} is the number of rows/columns in the covariance matrix \code{S},
            with the indices corresponding to the variables' positions in \code{S}.
            Variable indices above \eqn{n} represent latent variables.}
        \item{3. Arrow \emph{from}:}{the index of the variable at the tail of 
            a directional arrow, or at the other end of a bidirectional arrow.}
        \item{4. Parameter number:}{free parameters are numbered from 1 to \eqn{t},
            but do not necessarily appear in consecutive order. Fixed parameters are given
            the number 0. Equality contraints are specified by assigning two or more 
            parameters the same number.}
        \item{5. Value:}{start value for a free parameter, or value of a fixed parameter. If given
            as \code{NA}, the program will compute a start value, by a slight modification of the
            method described by McDonald and Hartmann (1992). \emph{Note:} In some circumstances,
            some start values are selected randomly; this might produce small differences in
            the parameter estimates when the program is rerun.}
        }
        
  The numeric \code{ram} matrix is normally generated automatically, not specified directly by the user.
        
  \code{sem} fits the model by calling the optimizer specified in the \code{optimizer} argument
  to minimize the objective function specified in the \code{objective} argument.
  If the optimization fails to converge, a warning message is printed.
  
  The RAM formulation of the general structural equation model is given by the basic equation
  \deqn{v = Av + u}
  where \eqn{v} and \eqn{u} are vectors of random variables (observed or unobserved), and
  the parameter matrix \eqn{A} contains regression coefficients, symbolized by single-headed arrows
  in a path diagram. Another parameter matrix, 
  \deqn{P = E(uu')}
  contains covariances among the elements of \eqn{u} (assuming that the elements of \eqn{u} have zero
  means). Usually \eqn{v} contains endogenous and exogenous observed and unobserved variables, but not
  error variables (see the examples below).
  
  The \code{startvalues} function may be called directly, but is usually called by \code{sem.default}.
}

\value{
  \code{sem} returns an object of class \code{c(}\emph{objective}\code{, "sem")}, where \emph{objective}
  	is the name of the objective function that was optimized (e.g., \code{"objectiveML"}), with the following elements:
  \item{var.names}{vector of variable names.}
  \item{ram}{RAM matrix, including any rows generated for covariances
    among fixed exogenous variables; column 5 includes computed start values.}
  \item{S}{observed covariance matrix.}
  \item{J}{RAM selection matrix, \eqn{J}, which picks out observed variables.}
  \item{n.fix}{number of fixed exogenous variables.}
  \item{n}{number of observed variables.}
  \item{N}{number of observations.}
  \item{m}{number of variables (observed plus unobserved).}
  \item{t}{number of free parameters.}
  \item{raw}{\code{TRUE} if the model is fit to a raw moment matrix, \code{FALSE} otherwise.}
  \item{data}{the observed-variable data matrix, or \code{NULL} if data are not supplied.}
  \item{semmod}{the \code{semmod} specification object for the model, if one was supplied; otherwise \code{NULL}.}
  \item{optimizer}{the optimizer function.}
  \item{objective}{the objective function.}
  \item{coeff}{estimates of free parameters.}
  \item{vcov}{estimated asymptotic covariance matrix of parameter estimates, based on a numeric Hessian,
  	if supplied by the optimizer; otherwise \code{NULL}.}
  \item{par.posn}{indices of free parameters.}
  \item{convergence}{\code{TRUE} or \code{FALSE}, depending upon whether the optimization apparently converged.}
  \item{iterations}{number of iterations performed.}
  \item{criterion}{value of the objective function at the minimum.}
  \item{C}{model-reproduced covariance matrix.}
  \item{A}{RAM \eqn{A} matrix.}
  \item{P}{RAM \eqn{P} matrix.}
  \item{adj.obj}{robust adjusted value of the objective function; \code{NULL} if
  	\code{robust} is \code{FALSE}.}
  \item{robust.vcov}{robust estimated coefficient covariance matrix; \code{NULL} if
  	\code{robust} is \code{FALSE}.}
}

\section{Warning}{
    A common error is to fail to specify variance or covariance terms in the model, which are denoted 
    by double-headed arrows, \code{<->}.
    
    In general, every observed or latent variable in the model should be associated 
    with a variance or error variance. 
    This may be
    a free parameter to estimate or a fixed constant (as in the case of a latent exogenous 
    variable for which you wish to fix the variance, e.g.,
    to 1). Again in general, there will be an \emph{error variance} associated with each 
    endogenous variable in the model (i.e., each variable
    to which at least one single-headed arrow points --- including observed indicators of latent variables), 
    and a \emph{variance} associated with each exogenous variable (i.e., each variable that
    appears only at the tail of single-headed arrows, never at the head).
    
    To my knowledge, the only \emph{apparent} exception to this rule is for observed variables that 
    are declared to be fixed exogenous variables.
    In this case, the program generates the necessary (fixed-constant) variances and covariances automatically.
    
    If there are missing variances, a warning message will be printed, and estimation will almost surely 
    fail in some manner. Missing
    variances might well indicate that there are missing covariances too, but it is not possible 
    to deduce this in a mechanical manner. The \code{\link{specifyModel}} funciton will by default supply
    error-variance parameters if these are missing.
}

\references{

	Fox, J. (2006)
	Structural equation modeling with the sem package in R.
	\emph{Structural Equation Modeling} 
	\bold{13}:465--486.
	
    Bollen, K. A. (1989) 
    \emph{Structural Equations With Latent Variables.} Wiley.
    
    Bollen, K. A. and Long, J. S. (eds.)
    \emph{Testing Structural Equation Models}, Sage.
    
    McArdle, J. J. and Epstein, D. (1987)
    Latent growth curves within developmental structural equation models.
    \emph{Child Development}
    \bold{58}, 110--133.
    
    McArdle, J. J. and McDonald, R. P. (1984)
    Some algebraic properties of the reticular action model.
    \emph{British Journal of Mathematical and Statistical Psychology}
    \bold{37}, 234--251.
    
    McDonald, R. P. and Hartmann, W. M. (1992)
    A procedure for obtaining initial values of parameters in
    the RAM model. \emph{Multivariate Behavioral Research}
    \bold{27}, 57--76.
    
    Raftery, A. E. (1993)
    Bayesian model selection in structural equation models.
    In Bollen, K. A. and Long, J. S. (eds.)
    \emph{Testing Structural Equation Models}, Sage.
    
    Raftery, A. E. (1995)
    Bayesian model selection in social research (with discussion).
    \emph{Sociological Methodology} \bold{25}, 111--196.
    
    Satorra, A. (2000)
	Scaled and adjusted restricted tests in multi-sample analysis of moment structures. 
	pp. 233--247 in Heijmans, R.D.H., Pollock, D.S.G. & Satorra, A. (eds.)
	\emph{Innovations in Multivariate Statistical Analysis. A Festschrift for Heinz Neudecker }, Kluwer. 
}

\author{John Fox \email{jfox@mcmaster.ca} and Jarrett Byrnes}

\seealso{\code{\link{rawMoments}}, \code{\link{startvalues}}, 
	\code{\link{objectiveML}}, \code{\link{objectiveGLS}},  
	\code{\link{optimizerNlm}}, \code{\link{optimizerOptim}}, \code{\link{optimizerNlminb}},
	\code{\link[stats]{nlm}}, \code{\link[stats]{optim}}, \code{\link[stats]{nlminb}}
}

\examples{

# Note: The first set of examples can't be run via example() because the default file
#  argument of specifyModel() and readMoments() requires that the model 
#  specification and covariances, correlations, or raw moments be entered
#  at the command prompt. The examples can be copied and run in the R console,
#  however. See ?specifyModel and ?readMoments for further information.
#  These examples are repeated below using file input to specifyModel() and
#  readMoments(). The second version of the examples may be executed through example().

    \dontrun{

# ------------- Duncan, Haller and Portes peer-influences model ----------------------
# A nonrecursive SEM with unobserved endogenous variables and fixed exogenous variables

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
            
# Fit the model using a symbolic ram specification

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
    
# an equivalent specification, allowing specifyModel() to generate
#  variance parameters for endogenous variables (and suppressing the
#  unnecessary NAs):
 
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
 
sem.dhp.1 <- sem(model.dhp, R.DHP, 329,
    fixed.x=c('RParAsp', 'RIQ', 'RSES', 'FSES', 'FIQ', 'FParAsp'))
summary(sem.dhp.1)

# Fit the model using a numerical ram specification (not recommended!)

ram.dhp <- matrix(c(
#               heads   to      from    param  start
                1,       1,     11,      0,     1,
                1,       2,     11,      1,     NA, # lam21
                1,       3,     12,      0,     1,
                1,       4,     12,      2,     NA, # lam42
                1,      11,      5,      3,     NA, # gam11
                1,      11,      6,      4,     NA, # gam12
                1,      11,      7,      5,     NA, # gam13
                1,      11,      8,      6,     NA, # gam14
                1,      12,      7,      7,     NA, # gam23
                1,      12,      8,      8,     NA, # gam24
                1,      12,      9,      9,     NA, # gam25
                1,      12,     10,     10,     NA, # gam26
                1,      11,     12,     11,     NA, # beta12
                1,      12,     11,     12,     NA, # beta21
                2,       1,      1,     13,     NA, # theta1
                2,       2,      2,     14,     NA, # theta2
                2,       3,      3,     15,     NA, # theta3
                2,       4,      4,     16,     NA, # theta4
                2,      11,     11,     17,     NA, # psi11
                2,      12,     12,     18,     NA, # psi22
                2,      11,     12,     19,     NA  # psi12
                ), ncol=5, byrow=TRUE)

params.dhp <- c('lam21', 'lam42', 'gam11', 'gam12', 'gam13', 'gam14',
                 'gam23',  'gam24',  'gam25',  'gam26',
                 'beta12', 'beta21', 'theta1', 'theta2', 'theta3', 'theta4',
                 'psi11', 'psi22', 'psi12')
                 
vars.dhp <- c('ROccAsp', 'REdAsp', 'FOccAsp', 'FEdAsp', 'RParAsp', 'RIQ',
                'RSES', 'FSES', 'FIQ', 'FParAsp', 'RGenAsp', 'FGenAsp')
                
sem.dhp.2 <- sem(ram.dhp, R.DHP, 329, param.names=params.dhp, var.names=vars.dhp, 
	fixed.x=5:10)
summary(sem.dhp.2)


# -------------------- Wheaton et al. alienation data ----------------------
    

S.wh <- readMoments(names=c('Anomia67','Powerless67','Anomia71',
                                    'Powerless71','Education','SEI'))
   11.834                                    
    6.947    9.364                            
    6.819    5.091   12.532                    
    4.783    5.028    7.495    9.986            
   -3.839   -3.889   -3.841   -3.625   9.610     
  -21.899  -18.831  -21.748  -18.775  35.522  450.288

# This is the model in the SAS manual for PROC CALIS: A Recursive SEM with
# latent endogenous and exogenous variables.
# Curiously, both factor loadings for two of the latent variables are fixed.

model.wh.1 <- specifyModel()
    Alienation67   ->  Anomia67,      NA,     1
    Alienation67   ->  Powerless67,   NA,     0.833
    Alienation71   ->  Anomia71,      NA,     1
    Alienation71   ->  Powerless71,   NA,     0.833 
    SES            ->  Education,     NA,     1     
    SES            ->  SEI,           lamb,   NA
    SES            ->  Alienation67,  gam1,   NA
    Alienation67   ->  Alienation71,  beta,   NA
    SES            ->  Alienation71,  gam2,   NA
    Anomia67       <-> Anomia67,      the1,   NA
    Anomia71       <-> Anomia71,      the1,   NA
    Powerless67    <-> Powerless67,   the2,   NA
    Powerless71    <-> Powerless71,   the2,   NA
    Education      <-> Education,     the3,   NA
    SEI            <-> SEI,           the4,   NA
    Anomia67       <-> Anomia71,      the5,   NA
    Powerless67    <-> Powerless71,   the5,   NA
    Alienation67   <-> Alienation67,  psi1,   NA
    Alienation71   <-> Alienation71,  psi2,   NA
    SES            <-> SES,           phi,    NA
                           
sem.wh.1 <- sem(model.wh.1, S.wh, 932)
summary(sem.wh.1)

# The same model, but treating one loading for each latent variable as free
# (and equal to each other).

model.wh.2 <- specifyModel()
    Alienation67   ->  Anomia67,      NA,        1
    Alienation67   ->  Powerless67,   lamby,    NA
    Alienation71   ->  Anomia71,      NA,        1
    Alienation71   ->  Powerless71,   lamby,    NA 
    SES            ->  Education,     NA,        1     
    SES            ->  SEI,           lambx,    NA
    SES            ->  Alienation67,  gam1,     NA
    Alienation67   ->  Alienation71,  beta,     NA
    SES            ->  Alienation71,  gam2,     NA
    Anomia67       <-> Anomia67,      the1,     NA
    Anomia71       <-> Anomia71,      the1,     NA
    Powerless67    <-> Powerless67,   the2,     NA
    Powerless71    <-> Powerless71,   the2,     NA
    Education      <-> Education,     the3,     NA
    SEI            <-> SEI,           the4,     NA
    Anomia67       <-> Anomia71,      the5,     NA
    Powerless67    <-> Powerless71,   the5,     NA
    Alienation67   <-> Alienation67,  psi1,     NA
    Alienation71   <-> Alienation71,  psi2,     NA
    SES            <-> SES,           phi,      NA 


sem.wh.2 <- sem(model.wh.2, S.wh, 932)
summary(sem.wh.2)

# Compare the two models by a likelihood-ratio test:

anova(sem.wh.1, sem.wh.2)


# ----------------------- Thurstone data ---------------------------------------
#  Second-order confirmatory factor analysis, from the SAS manual for PROC CALIS

R.thur <- readMoments(diag=FALSE, names=c('Sentences','Vocabulary',
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
            
model.thur <- specifyModel()
    F1 -> Sentences,                      lam11
    F1 -> Vocabulary,                     lam21
    F1 -> Sent.Completion,                lam31
    F2 -> First.Letters,                  lam41
    F2 -> 4.Letter.Words,                 lam52
    F2 -> Suffixes,                       lam62
    F3 -> Letter.Series,                  lam73
    F3 -> Pedigrees,                      lam83
    F3 -> Letter.Group,                   lam93
    F4 -> F1,                             gam1
    F4 -> F2,                             gam2
    F4 -> F3,                             gam3
    F1 <-> F1,                            NA,     1
    F2 <-> F2,                            NA,     1
    F3 <-> F3,                            NA,     1
    F4 <-> F4,                            NA,     1

sem.thur <- sem(model.thur, R.thur, 213)
summary(sem.thur)

#------------------------- Kerchoff/Kenney path analysis ---------------------
# An observed-variable recursive SEM from the LISREL manual

R.kerch <- readMoments(diag=FALSE, names=c('Intelligence','Siblings',
                        'FatherEd','FatherOcc','Grades','EducExp','OccupAsp'))
    -.100                                
     .277  -.152                          
     .250  -.108  .611                     
     .572  -.105  .294   .248               
     .489  -.213  .446   .410   .597         
     .335  -.153  .303   .331   .478   .651   
    
model.kerch <- specifyModel()
    Intelligence -> Grades,       gam51
    Siblings -> Grades,           gam52
    FatherEd -> Grades,           gam53
    FatherOcc -> Grades,          gam54
    Intelligence -> EducExp,      gam61
    Siblings -> EducExp,          gam62
    FatherEd -> EducExp,          gam63
    FatherOcc -> EducExp,         gam64
    Grades -> EducExp,            beta65
    Intelligence -> OccupAsp,     gam71
    Siblings -> OccupAsp,         gam72
    FatherEd -> OccupAsp,         gam73
    FatherOcc -> OccupAsp,        gam74
    Grades -> OccupAsp,           beta75
    EducExp -> OccupAsp,          beta76
                       
sem.kerch <- sem(model.kerch, R.kerch, 737, fixed.x=c('Intelligence','Siblings',
    'FatherEd','FatherOcc'))
summary(sem.kerch)


#------------------- McArdle/Epstein latent-growth-curve model -----------------
# This model, from McArdle and Epstein (1987, p.118), illustrates the use of a 
# raw moment matrix to fit a model with an intercept. (The example was suggested
# by Mike Stoolmiller.)

M.McArdle <- readMoments(names=c('WISC1', 'WISC2', 'WISC3', 'WISC4', 'UNIT'))
    365.661                                      
    503.175     719.905                           
    675.656     958.479    1303.392                
    890.680    1265.846    1712.475    2278.257     
     18.034      25.819      35.255      46.593     1.000
 
mod.McArdle <- specifyModel()
    C -> WISC1, NA, 6.07
    C -> WISC2, B2, NA
    C -> WISC3, B3, NA
    C -> WISC4, B4, NA
    UNIT -> C, Mc, NA
    C <-> C, Vc, NA,
    WISC1 <-> WISC1, Vd, NA
    WISC2 <-> WISC2, Vd, NA
    WISC3 <-> WISC3, Vd, NA
    WISC4 <-> WISC4, Vd, NA

sem.McArdle <- sem(mod.McArdle, M.McArdle, 204, fixed.x="UNIT", raw=TRUE)
summary(sem.McArdle)
    
#------------ Bollen industrialization and democracy example -----------------
# This model, from Bollen (1989, Ch. 8), illustrates the use in sem() of a
# case-by-variable data (see ?Bollen) set rather than a covariance or moment matrix

model.bollen <- specifyModel()
	Demo60 -> y1, NA, 1
	Demo60 -> y2, lam2, 
	Demo60 -> y3, lam3, 
	Demo60 -> y4, lam4, 
	Demo65 -> y5, NA, 1
	Demo65 -> y6, lam2, 
	Demo65 -> y7, lam3, 
	Demo65 -> y8, lam4, 
	Indust -> x1, NA, 1
	Indust -> x2, lam6, 
	Indust -> x3, lam7, 
	y1 <-> y5, theta15
	y2 <-> y4, theta24
	y2 <-> y6, theta26
	y3 <-> y7, theta37
	y4 <-> y8, theta48
	y6 <-> y8, theta68
	Indust -> Demo60, gamma11, 
	Indust -> Demo65, gamma21, 
	Demo60 -> Demo65, beta21, 
	Indust <-> Indust, phi
	
sem.bollen <- sem(model.bollen, data=Bollen)
summary(sem.bollen)
summary(sem.bollen, robust=TRUE) # robust SEs and tests
summary(sem.bollen, analytic.se=FALSE) # uses numeric rather than analytic Hessian

sem.bollen.gls <- sem(model.bollen, data=Bollen, objective=objectiveGLS) # GLS rather than ML estimator
summary(sem.bollen.gls)
	}

## ===============================================================================
	
# The following examples use file input and may be executed via example():

etc <- file.path(.path.package(package="sem")[1], "etc") # path to data and model files

# ------------- Duncan, Haller and Portes peer-influences model ----------------------
# A nonrecursive SEM with unobserved endogenous variables and fixed exogenous variables

(R.DHP <- readMoments(file=file.path(etc, "R-DHP.txt"),
				diag=FALSE, names=c("ROccAsp", "REdAsp", "FOccAsp", 
                "FEdAsp", "RParAsp", "RIQ", "RSES", "FSES", "FIQ", "FParAsp")))
(model.dhp <- specifyModel(file=file.path(etc, "model-DHP.txt")))
sem.dhp.1 <- sem(model.dhp, R.DHP, 329,
    fixed.x=c('RParAsp', 'RIQ', 'RSES', 'FSES', 'FIQ', 'FParAsp'))
summary(sem.dhp.1)


# -------------------- Wheaton et al. alienation data ----------------------

(S.wh <- readMoments(file=file.path(etc, "S-Wheaton.txt"),
					names=c('Anomia67','Powerless67','Anomia71',
                            'Powerless71','Education','SEI')))

# This is the model in the SAS manual for PROC CALIS: A Recursive SEM with
# latent endogenous and exogenous variables.
# Curiously, both factor loadings for two of the latent variables are fixed.

(model.wh.1 <- specifyModel(file=file.path(etc, "model-Wheaton-1.txt")))                    
sem.wh.1 <- sem(model.wh.1, S.wh, 932)
summary(sem.wh.1)

# The same model, but treating one loading for each latent variable as free
# (and equal to each other).

(model.wh.2 <- specifyModel(file=file.path(etc, "model-Wheaton-2.txt")))
sem.wh.2 <- sem(model.wh.2, S.wh, 932)
summary(sem.wh.2)

# Compare the two models by a likelihood-ratio test:

anova(sem.wh.1, sem.wh.2)


# ----------------------- Thurstone data ---------------------------------------

#  Second-order confirmatory factor analysis, from the SAS manual for PROC CALIS

(R.thur <- readMoments(file=file.path(etc, "R-Thurstone.txt"),
		diag=FALSE, names=c('Sentences','Vocabulary',
        'Sent.Completion','First.Letters','4.Letter.Words','Suffixes',
        'Letter.Series','Pedigrees', 'Letter.Group')))
(model.thur <- specifyModel(file=file.path(etc, "model-Thurstone.txt")))
sem.thur <- sem(model.thur, R.thur, 213)
summary(sem.thur)


#------------------------- Kerchoff/Kenney path analysis ---------------------

# An observed-variable recursive SEM from the LISREL manual

(R.kerch <- readMoments(file=file.path(etc, "R-Kerchoff.txt"),
					   diag=FALSE, names=c('Intelligence','Siblings',
                        'FatherEd','FatherOcc','Grades','EducExp','OccupAsp')))
(model.kerch <- specifyModel(file=file.path(etc, "model-Kerchoff.txt")))
sem.kerch <- sem(model.kerch, R.kerch, 737, fixed.x=c('Intelligence','Siblings',
    'FatherEd','FatherOcc'))
summary(sem.kerch)


#------------------- McArdle/Epstein latent-growth-curve model -----------------

# This model, from McArdle and Epstein (1987, p.118), illustrates the use of a 
# raw moment matrix to fit a model with an intercept. (The example was suggested
# by Mike Stoolmiller.)

(M.McArdle <- readMoments(file=file.path(etc, "M-McArdle.txt"),
						 names=c('WISC1', 'WISC2', 'WISC3', 'WISC4', 'UNIT')))
(mod.McArdle <- specifyModel(file=file.path(etc, "model-McArdle.txt")))
sem.McArdle <- sem(mod.McArdle, M.McArdle, 204, fixed.x="UNIT", raw=TRUE)
summary(sem.McArdle)


#------------ Bollen industrialization and democracy example -----------------

# This model, from Bollen (1989, Ch. 8), illustrates the use in sem() of a
# case-by-variable data set (see ?Bollen) rather than a covariance or moment matrix

(model.bollen <- specifyModel(file=file.path(etc, "model-Bollen.txt")))
sem.bollen <- sem(model.bollen, data=Bollen)
summary(sem.bollen)
summary(sem.bollen, robust=TRUE) # robust SEs and tests
summary(sem.bollen, analytic.se=FALSE) # uses numeric rather than analytic Hessian

sem.bollen.gls <- sem(model.bollen, data=Bollen, objective=objectiveGLS) # GLS rather than ML estimator
summary(sem.bollen.gls)
}   

\keyword{models}
