\name{sem} 

\alias{sem}
\alias{sem.mod}
\alias{sem.default}
\alias{startvalues}
\alias{print.sem}
\alias{summary.sem}
\alias{print.summary.sem}

\title{General Structural Equation Models}
\description{
  \code{sem} fits general structural equation models (with both observed and
  unobserved variables) by the method of maximum likelihood, assuming
  multinormal errors. Observed variables are also called \emph{indicators} or 
  \emph{manifest variables}; unobserved variables are also called \emph{factors}
  or \emph{latent variables}. Normally, the generic function (\code{sem}) would
  be used.
}

\usage{
sem(ram, ...)

sem.mod(ram, S, N, obs.variables=rownames(S), fixed.x=NULL, debug=F, ...)

sem.default(ram, S, N, param.names=paste("Param", 1:t, sep = ""), 
    var.names=paste("V", 1:m, sep = ""), fixed.x=NULL, 
    analytic.gradient=T, heywood=F, warn=F, control=list())
    
startvalues(S, ram)

print.sem(x, ...)

summary.sem(object, digits=5, ...)
}

\arguments{
  \item{ram}{RAM specification, which is a simple encoding of the path
    diagram for the model. The \code{ram} matrix may be given either in symbolic (character)
    form, invoking \code{sem.mod}, which calls \code{sem.default} after setting up the model,
    or (less conveniently) in numeric form, invoking \code{sem.default} directly
    (see \bold{Details} below).}
  \item{S}{covariance matrix among observed variables; may be input as a symmetric matrix,
    or as a lower- or upper-triangular matrix.}
  \item{N}{number of observations on which the covariance matrix is based.}
  \item{obs.variables}{names of observed variables, by default taken from the row names of
    the covariance matrix \code{S}.}
  \item{fixed.x}{names (if the \code{ram} matrix is given in symbolic form) or indices 
    (if it is in numeric form) of fixed exogenous variables. Specifying these obviates
    the necessity of having to fix the variances and covariances among these 
    variables (and produces correct degrees of freedom for the model chisquare).}
  \item{debug}{if \code{TRUE}, some information is printed to help you debug the symbolic
    model specification; for example, if a variable name is misspelled, \code{sem} will
    assume that the variable is a (new) latent variable. The default is \code{FALSE}.}
  \item{...}{arguments to be passed down to \code{sem.default} (or ignored).}
  \item{param.names}{names of the \eqn{t} free parameters, given in their numerical order;
    default names are \code{Param1}, \ldots, \code{Param}\emph{t}. Note: Should not be
    specified when the \code{ram} matrix is given in symbolic form.}
  \item{var.names}{names of the \eqn{m} entries of the \eqn{v} vector
    (typically the observed and latent variables --- see below), given in their
    numerical order; default names are \code{Var1}, \ldots, \code{Var}\emph{m}.
    Note: Should not be specified when the \code{ram} matrix is given in symbolic form.}
  \item{analytic.gradient}{if \code{TRUE} (the default), then analytic first derivatives are
    used in the maximization of the likelihood; otherwise numeric derivatives are used.}
  \item{heywood}{if \code{TRUE}, variances are constrained to be non-negative in the 
    initial fitting process. The default is \code{FALSE}.}
  \item{warn}{if \code{TRUE}, warnings produced by the optimization functions will be printed.
    This should generally not be necessary, since \code{sem} prints its own warnings, and saves
    information about convergence. The default is \code{FALSE}.}
  \item{control}{a list with up to three elements:
    \code{optim.control}, a list of control parameters to be passed to 
            the \code{optim} function via its \code{control} argument;
    \code{optim.method}, the optimization method to be employed by \code{optim}; 
            if \code{heywood} is \code{TRUE}, then this defaults to
            \code{"L-BFGS-B"}, otherwise to \code{"BFGS"}.
    \code{nlm.iterlim}, the maximum number of iterations for the 
            \code{nlm} function, to be passed to it via its \code{iterlim}
            argument.}
  \item{object, x}{an object of class \code{sem} returned by the \code{sem} function.}
  \item{digits}{number of digits for printed output.}
} 

\details{
  The model is set up using RAM (`reticular action model' -- don't ask!) 
  notation -- a simple format
  for specifying general structural equation models by coding the
  `arrows' in the path diagram for the model (see, e.g., McArdle and McDonald, 1984).
  
  The variables in the \eqn{v} vector in the model (typically, the observed and
  unobserved variables, but not error variables) are numbered from 1 to \eqn{m}.  
  the RAM matrix contains one row for each (free or constrained) parameter of the model, and
  may be specified either in symbolic (character) format or in numeric format.
  
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
        
  \code{sem} fits the model by calling the \code{\link[base]{optim}}
  and \code{\link[base]{nlm}} optimizers
  sequentially to minimize the negative log-likelihood for the model. 
  If \code{nlm} fails to converge using an analytic gradient, then a numeric gradient is tried
  (and vice-versa if the numeric gradient is tried first). Under these circumstances, a warning
  message is printed.
  
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
  \code{sem} returns an object of class \code{sem}, with the following elements:
  \item{ram}{RAM matrix, including any rows generated for covariances
    among fixed exogenous variables; column 5 includes computed start values.}
  \item{coeff}{estimates of free parameters.}
  \item{criterion}{fitting criterion --- negative log-liklihood divided by \eqn{N - 1}.}
  \item{cov}{estimated asymptotic covariance matrix of parameter estimates.}
  \item{S}{observed covariance matrix.}
  \item{J}{RAM selection matrix, \eqn{J}, which picks out observed variables.}
  \item{C}{model-reproduced covariance matrix.}
  \item{A}{RAM \eqn{A} matrix.}
  \item{P}{RAM \eqn{P} matrix.}
  \item{n.fix}{number of fixed exogenous variables.}
  \item{n}{number of observed variables.}
  \item{N}{number of observations.}
  \item{m}{number of variables (observed plus unobserved).}
  \item{t}{number of free parameters.}
  \item{par.posn}{indices of free parameters.}
  \item{var.names}{vector of variable names.}
  \item{observed}{indices of observed variables.}
  \item{convergence.1}{convergence code returned by \code{optim} (a code > 0
    indicates a problem).}
  \item{message.1}{message returned by \code{optim} in the first stage of
    optimization.}
  \item{coef.1}{parameter estimates returned by \code{optim}.}
  \item{convergence.2}{convergence code returned by \code{nlm} in the
    second stage of optimization (a code > 2 indicates a problem).}
  \item{coef.2}{parameter estimates returned by \code{nlm} in
    second stage of optimization. Present only if second stage failed.}
  \item{convergence.3}{convergence code returned by \code{nlm} in the
    third stage of optimization. Present only if second stage failed.}
}

\references{
    Bollen, K. A. (1989) 
    \emph{Structural Equations With Latent Variables.} Wiley.
    
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
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[base]{optim}}, \code{\link[base]{nlm}}}

\examples{

# ------------- Duncan, Haller and Portes peer-influences model ----------------------
# A nonrecursive SEM with unobserved endogenous variables and fixed exogenous variables

R.DHP <- matrix(c(      # lower triangle of correlation matrix
            1,      0,      0,      0,      0,      0,      0,      0,      0,      0,            
            .6247,   1,     0,      0,      0,      0,      0,      0,      0,      0,            
            .3269,  .3669,  1,      0,      0,      0,      0,      0,      0,      0,            
            .4216,  .3275,  .6404,  1,      0,      0,      0,      0,      0,      0,
            .2137,  .2742,  .1124,  .0839,  1,      0,      0,      0,      0,      0,
            .4105,  .4043,  .2903,  .2598,  .1839,  1,      0,      0,      0,      0,
            .3240,  .4047,  .3054,  .2786,  .0489,  .2220,  1,      0,      0,      0,
            .2930,  .2407,  .4105,  .3607,  .0186,  .1861,  .2707,  1,      0,      0,
            .2995,  .2863,  .5191,  .5007,  .0782,  .3355,  .2302,  .2950,  1,      0,
            .0760,  .0702,  .2784,  .1988,  .1147,  .1021,  .0931, -.0438,  .2087,  1
            ), ncol=10, byrow=T)

# Fit the model using a symbolic ram specification

model.dhp <- matrix(c(
                      'RParAsp  -> RGenAsp', 'gam11',  NA,
                      'RIQ      -> RGenAsp', 'gam12',  NA,
                      'RSES     -> RGenAsp', 'gam13',  NA,
                      'FSES     -> RGenAsp', 'gam14',  NA,
                      'RSES     -> FGenAsp', 'gam23',  NA,
                      'FSES     -> FGenAsp', 'gam24',  NA,
                      'FIQ      -> FGenAsp', 'gam25',  NA,
                      'FParAsp  -> FGenAsp', 'gam26',  NA,
                      'FGenAsp  -> RGenAsp', 'beta12', NA,
                      'RGenAsp  -> FGenAsp', 'beta21', NA,
                      'RGenAsp  -> ROccAsp',  NA,       1,
                      'RGenAsp  -> REdAsp',  'lam21',  NA,
                      'FGenAsp  -> FOccAsp',  NA,       1,
                      'FGenAsp  -> FEdAsp',  'lam42',  NA,
                      'RGenAsp <-> RGenAsp', 'ps11',   NA,
                      'FGenAsp <-> FGenAsp', 'ps22',   NA,
                      'RGenAsp <-> FGenAsp', 'ps12',   NA,
                      'ROccAsp <-> ROccAsp', 'theta1', NA,
                      'REdAsp  <-> REdAsp',  'theta2', NA,
                      'FOccAsp <-> FOccAsp', 'theta3', NA,
                      'FEdAsp  <-> FEdAsp',  'theta4', NA),
                    ncol=3, byrow=T)

obs.vars.dhp <- c('ROccAsp', 'REdAsp', 'FOccAsp', 'FEdAsp', 'RParAsp', 'RIQ',
                'RSES', 'FSES', 'FIQ', 'FParAsp')

sem.dhp.1 <- sem(model.dhp, R.DHP, 329, obs.vars.dhp, 
    fixed.x=c('RParAsp', 'RIQ', 'RSES', 'FSES', 'FIQ', 'FParAsp'))

summary(sem.dhp.1)

##  Model Chisquare =  26.697   Df =  15 Pr(>Chisq) = 0.031302
##  Goodness-of-fit index =  0.98439
##  Adjusted goodness-of-fit index =  0.94275
##  BIC =  -94.782 
##  
##  Normalized Residuals
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -8.01e-01 -1.18e-01  5.02e-16 -1.20e-02  3.98e-02  1.57e+00 
##  
##  Parameter Estimates
##          Estimate Std Error  z value   Pr(>|z|)                     
##  gam11   0.161222  0.038487  4.18896 1.4012e-05 RGenAsp <--- RParAsp
##  gam12   0.249649  0.044581  5.59995 1.0720e-08     RGenAsp <--- RIQ
##  gam13   0.218404  0.043476  5.02355 2.5362e-07    RGenAsp <--- RSES
##  gam14   0.071841  0.050335  1.42724 7.6755e-02    RGenAsp <--- FSES
##  gam23   0.061885  0.051738  1.19612 1.1583e-01    FGenAsp <--- RSES
##  gam24   0.228867  0.044495  5.14368 1.3471e-07    FGenAsp <--- FSES
##  gam25   0.349036  0.044551  7.83455 2.3315e-15     FGenAsp <--- FIQ
##  gam26   0.159535  0.040129  3.97552 3.5113e-05 FGenAsp <--- FParAsp
##  beta12  0.184227  0.096207  1.91490 2.7753e-02 RGenAsp <--- FGenAsp
##  beta21  0.235487  0.119744  1.96659 2.4616e-02 FGenAsp <--- RGenAsp
##  lam21   1.062682  0.091969 11.55481 0.0000e+00  REdAsp <--- RGenAsp
##  lam42   0.929717  0.071151 13.06687 0.0000e+00  FEdAsp <--- FGenAsp
##  ps11    0.280987  0.046311  6.06740 6.5001e-10 RGenAsp <--> RGenAsp
##  ps22    0.263838  0.044902  5.87587 2.1032e-09 FGenAsp <--> FGenAsp
##  ps12   -0.022611  0.051650 -0.43778 3.3077e-01 FGenAsp <--> RGenAsp
##  theta1  0.412147  0.052211  7.89389 1.4433e-15 ROccAsp <--> ROccAsp
##  theta2  0.336144  0.053323  6.30390 1.4512e-10   REdAsp <--> REdAsp
##  theta3  0.311191  0.046665  6.66866 1.2907e-11 FOccAsp <--> FOccAsp
##  theta4  0.404607  0.046733  8.65782 0.0000e+00   FEdAsp <--> FEdAsp


\dontrun{

# Fit the model using a numerical ram specification

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
                ), ncol=5, byrow=T)

params.dhp <- c('lam21', 'lam42', 'gam11', 'gam12', 'gam13', 'gam14',
                 'gam23',  'gam24',  'gam25',  'gam26',
                 'beta12', 'beta21', 'theta1', 'theta2', 'theta3', 'theta4',
                 'psi11', 'psi22', 'psi12')
                 
vars.dhp <- c('ROccAsp', 'REdAsp', 'FOccAsp', 'FEdAsp', 'RParAsp', 'RIQ',
                'RSES', 'FSES', 'FIQ', 'FParAsp', 'RGenAsp', 'FGenAsp')
                
sem.dhp.2 <- sem(ram.dhp, R.DHP, 329, params.dhp, vars.dhp, fixed.x=5:10)

summary(sem.dhp.2)
    
##  Model Chisquare =  26.697   Df =  15 Pr(>Chisq) = 0.031302
##  Goodness-of-fit index =  0.98439
##  Adjusted goodness-of-fit index =  0.94275
##  BIC =  -94.782 
##  
##  Normalized Residuals
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -8.01e-01 -1.18e-01  7.53e-16 -1.20e-02  3.98e-02  1.57e+00 
##  
##  Parameter Estimates
##          Estimate Std Error  z value   Pr(>|z|)                     
##  lam21   1.062682  0.091969 11.55481 0.0000e+00  REdAsp <--- RGenAsp
##  lam42   0.929717  0.071151 13.06687 0.0000e+00  FEdAsp <--- FGenAsp
##  gam11   0.161222  0.038487  4.18896 1.4012e-05 RGenAsp <--- RParAsp
##  gam12   0.249649  0.044581  5.59995 1.0720e-08     RGenAsp <--- RIQ
##  gam13   0.218404  0.043476  5.02355 2.5362e-07    RGenAsp <--- RSES
##  gam14   0.071841  0.050335  1.42724 7.6755e-02    RGenAsp <--- FSES
##  gam23   0.061885  0.051738  1.19612 1.1583e-01    FGenAsp <--- RSES
##  gam24   0.228867  0.044495  5.14368 1.3471e-07    FGenAsp <--- FSES
##  gam25   0.349036  0.044551  7.83455 2.3315e-15     FGenAsp <--- FIQ
##  gam26   0.159535  0.040129  3.97552 3.5113e-05 FGenAsp <--- FParAsp
##  beta12  0.184227  0.096207  1.91490 2.7753e-02 RGenAsp <--- FGenAsp
##  beta21  0.235487  0.119744  1.96658 2.4616e-02 FGenAsp <--- RGenAsp
##  theta1  0.412147  0.052211  7.89389 1.4433e-15 ROccAsp <--> ROccAsp
##  theta2  0.336144  0.053323  6.30390 1.4512e-10   REdAsp <--> REdAsp
##  theta3  0.311191  0.046665  6.66866 1.2907e-11 FOccAsp <--> FOccAsp
##  theta4  0.404607  0.046733  8.65782 0.0000e+00   FEdAsp <--> FEdAsp
##  psi11   0.280987  0.046311  6.06739 6.5001e-10 RGenAsp <--> RGenAsp
##  psi22   0.263838  0.044902  5.87587 2.1032e-09 FGenAsp <--> FGenAsp
##  psi12  -0.022611  0.051650 -0.43778 3.3077e-01 RGenAsp <--> FGenAsp



# -------------------- Wheaton et al. alienation data ----------------------
    

S.wh <- matrix(c(
   11.834,     0,        0,        0,       0,        0,
    6.947,    9.364,     0,        0,       0,        0,
    6.819,    5.091,   12.532,     0,       0,        0,
    4.783,    5.028,    7.495,    9.986,    0,        0,
   -3.839,   -3.889,   -3.841,   -3.625,   9.610,     0,
  -21.899,  -18.831,  -21.748,  -18.775,  35.522,  450.288), 
  6, 6)
  
# This is the model in the SAS manual for PROC CALIS: A Recursive SEM with
# latent endogenous and exogenous variables.
# Curiously, both factor loadings for two of the latent variables are fixed.

model.wh.1 <- matrix(c(
                    'Alienation67   ->  Anomia67',      NA,     1,
                    'Alienation67   ->  Powerless67',   NA,     0.833,
                    'Alienation71   ->  Anomia71',      NA,     1,
                    'Alienation71   ->  Powerless71',   NA,     0.833, 
                    'SES            ->  Education',     NA,     1,     
                    'SES            ->  SEI',           'lamb', NA,
                    'SES            ->  Alienation67',  'gam1', NA,
                    'Alienation67   ->  Alienation71',  'beta', NA,
                    'SES            ->  Alienation71',  'gam2', NA,
                    'Anomia67       <-> Anomia67',      'the1', NA,
                    'Anomia71       <-> Anomia71',      'the1', NA,
                    'Powerless67    <-> Powerless67',   'the2', NA,
                    'Powerless71    <-> Powerless71',   'the2', NA,
                    'Education      <-> Education',     'the3', NA,
                    'SEI            <-> SEI',           'the4', NA,
                    'Anomia67       <-> Anomia71',      'the5', NA,
                    'Powerless67    <-> Powerless71',   'the5', NA,
                    'Alienation67   <-> Alienation67',  'psi1', NA,
                    'Alienation71   <-> Alienation71',  'psi2', NA,
                    'SES            <-> SES',           'phi',  NA), 
                    ncol=3, byrow=T)
                        
obs.vars.wh <- c('Anomia67','Powerless67','Anomia71','Powerless71','Education','SEI')

sem.wh.1 <- sem(model.wh.1, S.wh, 932, obs.vars.wh)
##  Warning messages: 
##  1: second optimization DID NOT converge in: 
##      sem.default(ram = ram, S = S, N = N, param.names = pars, var.names = vars,  
##  2: final optimization DID converge in: 
##      sem.default(ram = ram, S = S, N = N, param.names = pars, var.names = vars,  

# Note failure of second optimization, but success of third

summary(sem.wh.1)

##  Model Chisquare =  13.485   Df =  9 Pr(>Chisq) = 0.14186
##  Goodness-of-fit index =  0.99527
##  Adjusted goodness-of-fit index =  0.98896
##  BIC =  -64.177 
##  
##  Normalized Residuals
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -1.26e+00 -1.31e-01  9.69e-05 -2.87e-02  1.14e-01  8.75e-01 
##  
##  Parameter Estimates
##      Estimate Std Error  z value   Pr(>|z|)                               
##  lamb   5.36885  0.433991  12.3709 0.0000e+00                   SEI <--- SES
##  gam1  -0.62994  0.056128 -11.2232 0.0000e+00          Alienation67 <--- SES
##  beta   0.59312  0.046820  12.6680 0.0000e+00 Alienation71 <--- Alienation67
##  gam2  -0.24086  0.055202  -4.3632 6.4085e-06          Alienation71 <--- SES
##  the1   3.60790  0.200591  17.9863 0.0000e+00         Anomia67 <--> Anomia67
##  the2   3.59496  0.165236  21.7566 0.0000e+00   Powerless67 <--> Powerless67
##  the3   2.99366  0.498977   5.9996 9.8905e-10       Education <--> Education
##  the4 259.57742 18.321464  14.1679 0.0000e+00                   SEI <--> SEI
##  the5   0.90580  0.121712   7.4422 4.9516e-14         Anomia71 <--> Anomia67
##  psi1   5.67051  0.422907  13.4084 0.0000e+00 Alienation67 <--> Alienation67
##  psi2   4.51485  0.334999  13.4772 0.0000e+00 Alienation71 <--> Alienation71
##  phi    6.61627  0.639511  10.3458 0.0000e+00                   SES <--> SES


# The same model, but treating one loading for each latent variable as free.

model.wh.2 <- matrix(c(
                    'Alienation67   ->  Anomia67',      NA,         1,
                    'Alienation67   ->  Powerless67',   'lamby',    NA,
                    'Alienation71   ->  Anomia71',      NA,         1,
                    'Alienation71   ->  Powerless71',   'lamby',    NA, 
                    'SES            ->  Education',     NA,         1,     
                    'SES            ->  SEI',           'lambx',    NA,
                    'SES            ->  Alienation67',  'gam1',     NA,
                    'Alienation67   ->  Alienation71',  'beta',     NA,
                    'SES            ->  Alienation71',  'gam2',     NA,
                    'Anomia67       <-> Anomia67',      'the1',     NA,
                    'Anomia71       <-> Anomia71',      'the1',     NA,
                    'Powerless67    <-> Powerless67',   'the2',     NA,
                    'Powerless71    <-> Powerless71',   'the2',     NA,
                    'Education      <-> Education',     'the3',     NA,
                    'SEI            <-> SEI',           'the4',     NA,
                    'Anomia67       <-> Anomia71',      'the5',     NA,
                    'Powerless67    <-> Powerless71',   'the5',     NA,
                    'Alienation67   <-> Alienation67',  'psi1',     NA,
                    'Alienation71   <-> Alienation71',  'psi2',     NA,
                    'SES            <-> SES',           'phi',      NA), 
                    ncol=3, byrow=T)                        

sem.wh.2 <- sem(model.wh.2, S.wh, 932, obs.vars.wh)
##  Warning messages: 
##  1: second optimization DID NOT converge in: 
##      sem.default(ram = ram, S = S, N = N, param.names = pars, var.names = vars,  
##  2: final optimization DID converge in: 
##      sem.default(ram = ram, S = S, N = N, param.names = pars, var.names = vars,  

summary(sem.wh.2)

##  Model Chisquare =  12.673   Df =  8 Pr(>Chisq) = 0.12360
##  Goodness-of-fit index =  0.99553
##  Adjusted goodness-of-fit index =  0.98828
##  BIC =  -56.36 
##  
##  Normalized Residuals
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -0.997000 -0.139000  0.000142 -0.028500  0.100000  0.759000 
##  
##  Parameter Estimates
##      Estimate Std Error  z value   Pr(>|z|)                               
##  lamby   0.86262  0.033383  25.8399 0.0000e+00  Powerless67 <--- Alienation67
##  lambx   5.35329  0.432588  12.3750 0.0000e+00                   SEI <--- SES
##  gam1   -0.62132  0.056144 -11.0666 0.0000e+00          Alienation67 <--- SES
##  beta    0.59424  0.047041  12.6324 0.0000e+00 Alienation71 <--- Alienation67
##  gam2   -0.23584  0.054688  -4.3126 8.0685e-06          Alienation71 <--- SES
##  the1    3.74502  0.249826  14.9905 0.0000e+00         Anomia67 <--> Anomia67
##  the2    3.49376  0.200757  17.4030 0.0000e+00   Powerless67 <--> Powerless67
##  the3    2.97442  0.499623   5.9533 1.3137e-09       Education <--> Education
##  the4  260.12964 18.297803  14.2164 0.0000e+00                   SEI <--> SEI
##  the5    0.90378  0.121818   7.4191 5.8953e-14         Anomia71 <--> Anomia67
##  psi1    5.47354  0.464054  11.7951 0.0000e+00 Alienation67 <--> Alienation67
##  psi2    4.36409  0.362720  12.0315 0.0000e+00 Alienation71 <--> Alienation71
##  phi     6.63548  0.640408  10.3613 0.0000e+00                   SES <--> SES


# ----------------------- Thurstone data ---------------------------------------
#  Second-order confirmatory factor analysis, from the SAS manual for PROC CALIS

R.thur <- matrix(c(
        1.,       0,      0,      0,      0,      0,      0,      0,      0, 
         .828,   1.,      0,      0,      0,      0,      0,      0,      0, 
         .776,   .779,   1.,      0,      0,      0,      0,      0,      0, 
         .439,   .493,    .460,   1.,     0,      0,      0,      0,      0, 
         .432,   .464,    .425,   .674,   1.,     0,      0,      0,      0, 
         .447,   .489,    .443,   .590,    .541,  1.,     0,      0,      0, 
         .447,   .432,    .401,   .381,    .402,   .288,  1.,     0,      0, 
         .541,   .537,    .534,   .350,    .367,   .320,   .555,  1.,     0, 
         .380,   .358,    .359,   .424,    .446,   .325,   .598,   .452,  1.
            ), ncol=9, byrow=T)

model.thur <- matrix(c(
                        'F1 -> Sentences',                      'lam11', NA,
                        'F1 -> Vocabulary',                     'lam21', NA,
                        'F1 -> Sent.Completion',                'lam31', NA,
                        'F2 -> First.Letters',                  'lam41', NA,
                        'F2 -> 4.Letter.Words',                  'lam52', NA,
                        'F2 -> Suffixes',                       'lam62', NA,
                        'F3 -> Letter.Series',                  'lam73', NA,
                        'F3 -> Pedigrees',                      'lam83', NA,
                        'F3 -> Letter.Group',                   'lam93', NA,
                        'F4 -> F1',                             'gam1',  NA,
                        'F4 -> F2',                             'gam2',  NA,
                        'F4 -> F3',                             'gam3',  NA, 
                        'Sentences <-> Sentences',              'th1',   NA,
                        'Vocabulary <-> Vocabulary',            'th2',   NA,
                        'Sent.Completion <-> Sent.Completion',  'th3',   NA,
                        'First.Letters <-> First.Letters',      'th4',   NA,
                        '4.Letter.Words <-> 4.Letter.Words',    'th5',   NA,
                        'Suffixes <-> Suffixes',                'th6',   NA,
                        'Letter.Series <-> Letter.Series',      'th7',   NA,
                        'Pedigrees <-> Pedigrees',              'th8',   NA,
                        'Letter.Group <-> Letter.Group',        'th9',   NA,
                        'F1 <-> F1',                            NA,      1,
                        'F2 <-> F2',                            NA,      1,
                        'F3 <-> F3',                            NA,      1,
                        'F4 <-> F4',                            NA,      1),
                        ncol=3, byrow=T)
                        
obs.vars.thur <- c('Sentences','Vocabulary','Sent.Completion','First.Letters',
                '4.Letter.Words','Suffixes','Letter.Series','Pedigrees',
                'Letter.Group')

sem.thur <- sem(model.thur, R.thur, 213, obs.vars.thur)

summary(sem.thur)

##  Model Chisquare =  38.196   Df =  24 Pr(>Chisq) = 0.033101
##  Goodness-of-fit index =  0.95957
##  Adjusted goodness-of-fit index =  0.9242
##  BIC =  -143.21 
##  
##  Normalized Residuals
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -9.75e-01 -4.18e-01 -6.45e-05  4.01e-02  9.40e-02  1.63e+00 
##  
##  Parameter Estimates
##      Estimate Std Error z value   Pr(>|z|)                                     
##  lam11  0.51513  0.064963  7.9296 1.1102e-15                    Sentences <--- F1
##  lam21  0.52032  0.065160  7.9852 6.6613e-16                   Vocabulary <--- F1
##  lam31  0.48744  0.062420  7.8090 2.8866e-15              Sent.Completion <--- F1
##  lam41  0.52110  0.063140  8.2532 1.1102e-16                First.Letters <--- F2
##  lam52  0.49705  0.059676  8.3292 0.0000e+00               4.Letter.Words <--- F2
##  lam62  0.43805  0.056481  7.7557 4.4409e-15                     Suffixes <--- F2
##  lam73  0.45243  0.071371  6.3391 1.1553e-10                Letter.Series <--- F3
##  lam83  0.41729  0.061036  6.8367 4.0515e-12                    Pedigrees <--- F3
##  lam93  0.40763  0.064524  6.3175 1.3295e-10                 Letter.Group <--- F3
##  gam1   1.44380  0.264162  5.4656 2.3071e-08                           F1 <--- F4
##  gam2   1.25391  0.216623  5.7884 3.5521e-09                           F2 <--- F4
##  gam3   1.40657  0.279336  5.0354 2.3841e-07                           F3 <--- F4
##  th1    0.18150  0.028401  6.3907 8.2580e-11             Sentences <--> Sentences
##  th2    0.16493  0.027797  5.9334 1.4840e-09           Vocabulary <--> Vocabulary
##  th3    0.26713  0.033469  7.9816 7.7716e-16 Sent.Completion <--> Sent.Completion
##  th4    0.30150  0.050686  5.9485 1.3531e-09     First.Letters <--> First.Letters
##  th5    0.36450  0.052358  6.9617 1.6803e-12   4.Letter.Words <--> 4.Letter.Words
##  th6    0.50641  0.059962  8.4455 0.0000e+00               Suffixes <--> Suffixes
##  th7    0.39033  0.061598  6.3367 1.1735e-10     Letter.Series <--> Letter.Series
##  th8    0.48137  0.065388  7.3618 9.0705e-14             Pedigrees <--> Pedigrees
##  th9    0.50510  0.065227  7.7437 4.7740e-15       Letter.Group <--> Letter.Group


#------------------------- Kerchoff/Kenney path analysis ---------------------
# An observed-variable recursive SEM from the LISREL manual

R.kerch <- matrix(c(
    1,      0,      0,      0,      0,      0,      0,
    -.100,  1,      0,      0,      0,      0,      0,
     .277,  -.152,  1,      0,      0,      0,      0,
     .250,  -.108,  .611,   1,      0,      0,      0,
     .572,  -.105,  .294,   .248,   1,      0,      0,
     .489,  -.213,  .446,   .410,   .597,   1,      0,
     .335,  -.153,  .303,   .331,   .478,   .651,   1),
     ncol=7, byrow=T)

rownames(R.kerch) <- colnames(R.kerch) <- c('Intelligence','Siblings',
    'FatherEd','FatherOcc','Grades','EducExp','OccupAsp')

    
model.kerch <- matrix(c(
                        'Intelligence -> Grades',       'gam51',    NA,
                        'Siblings -> Grades',           'gam52',    NA,
                        'FatherEd -> Grades',           'gam53',    NA,
                        'FatherOcc -> Grades',          'gam54',    NA,
                        'Intelligence -> EducExp',      'gam61',    NA,
                        'Siblings -> EducExp',          'gam62',    NA,
                        'FatherEd -> EducExp',          'gam63',    NA,
                        'FatherOcc -> EducExp',         'gam64',    NA,
                        'Grades -> EducExp',            'beta65',   NA,
                        'Intelligence -> OccupAsp',     'gam71',    NA,
                        'Siblings -> OccupAsp',         'gam72',    NA,
                        'FatherEd -> OccupAsp',         'gam73',    NA,
                        'FatherOcc -> OccupAsp',        'gam74',    NA,
                        'Grades -> OccupAsp',           'beta75',   NA,
                        'EducExp -> OccupAsp',          'beta76',   NA,
                        'Grades <-> Grades',            'psi5',     NA,
                        'EducExp <-> EducExp',          'psi6',     NA,
                        'OccupAsp <-> OccupAsp',        'psi7',     NA),
                       ncol=3, byrow=T)
                       
sem.kerch <- sem(model.kerch, R.kerch, 737, fixed.x=c('Intelligence','Siblings',
    'FatherEd','FatherOcc'))
    
summary(sem.kerch)

##  Model Chisquare =  6.537e-13   Df =  0 Pr(>Chisq) = NA
##  Goodness-of-fit index =  1
##  Adjusted goodness-of-fit index =  NA
##  BIC =  NA 
##  
##  Normalized Residuals
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -1.70e-14 -3.75e-16  7.45e-16  3.53e-15  4.47e-15  3.62e-14 
##  
##  Parameter Estimates
##          Estimate Std Error  z value   Pr(>|z|)                           
##  gam51   0.525902  0.031182 16.86530 0.0000e+00   Grades <--- Intelligence
##  gam52  -0.029942  0.030149 -0.99314 1.6032e-01       Grades <--- Siblings
##  gam53   0.118966  0.038259  3.10951 9.3699e-04       Grades <--- FatherEd
##  gam54   0.040603  0.037785  1.07456 1.4129e-01      Grades <--- FatherOcc
##  gam61   0.160270  0.032710  4.89979 4.7970e-07  EducExp <--- Intelligence
##  gam62  -0.111779  0.026876 -4.15899 1.5983e-05      EducExp <--- Siblings
##  gam63   0.172719  0.034306  5.03461 2.3941e-07      EducExp <--- FatherEd
##  gam64   0.151852  0.033688  4.50758 3.2785e-06     EducExp <--- FatherOcc
##  beta65  0.405150  0.032838 12.33799 0.0000e+00        EducExp <--- Grades
##  gam71  -0.039405  0.034500 -1.14215 1.2670e-01 OccupAsp <--- Intelligence
##  gam72  -0.018825  0.028222 -0.66700 2.5238e-01     OccupAsp <--- Siblings
##  gam73  -0.041333  0.036216 -1.14126 1.2688e-01     OccupAsp <--- FatherEd
##  gam74   0.099577  0.035446  2.80924 2.4829e-03    OccupAsp <--- FatherOcc
##  beta75  0.157912  0.037443  4.21738 1.2358e-05       OccupAsp <--- Grades
##  beta76  0.549593  0.038260 14.36486 0.0000e+00      OccupAsp <--- EducExp
##  psi5    0.650995  0.033946 19.17743 0.0000e+00         Grades <--> Grades
##  psi6    0.516652  0.026943 19.17590 0.0000e+00       EducExp <--> EducExp
##  psi7    0.556617  0.029026 19.17644 0.0000e+00     OccupAsp <--> OccupAsp

    }
}

\keyword{models}
