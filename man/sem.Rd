\name{sem}
\alias{sem}
\alias{start.values}
\alias{print.sem}
\alias{summary.sem}

\title{General Structural Equation Models}
\description{
  \code{sem} fits general structural equation models (with both observed and
  unobserved variables) by the method of maximum likelihood, assuming
  multinormal errors. Observed variables are also called \emph{indicators} or 
  \emph{manifest variables}; unobserved variables are also called \emph{factors}
  or \emph{latent variables}.
}
\usage{
sem(S, ram, N, param.names=paste("Param", 1:t, sep = ""), 
    var.names=paste("V", 1:m, sep=""), observed=1:n, 
    fixed.x=NULL, heywood=T, control=list())
    
start.values(S, ram, observed)

print.sem(x)

summary.sem(object, digits=5)
}

\arguments{
  \item{S}{covariance matrix among observed variables.}
  \item{ram}{RAM specification, which is a simple encoding of the path
    diagram for the model (see details below).}
  \item{N}{number of observations on which the covariance matrix is based.}
  \item{param.names}{names of the \eqn{t} free parameters, given in their numerical order;
    default names are \code{Param1}, \ldots, \code{Param}\emph{t}.}
  \item{var.names}{names of the \eqn{m} entries of the \eqn{v} vector
    (typically the observed and latent variables --- see below), given in their
    numerical order; default names are \code{Var1}, \ldots, \code{Var}\emph{m}.}
  \item{observed}{indices of the \eqn{n} observed variables; default is 1 to \eqn{n}.}
  \item{fixed.x}{indices of fixed exogenous variables. Specifying these obviates
    the necessity of having to fix the variances and covariances among these 
    variables (and produces correct degrees of freedom for the model chisquare).}
  \item{heywood}{if \code{TRUE} (the default) variances are constrained to be
    non-negative in the initial fitting process.}
  \item{control}{a list with up to two elements:
    \code{optim.control}, a list of control parameters to be passed to 
            the \code{optim} function via its \code{control} argument;
    \code{nlm.iterlim}, the maximum number of iterations for the 
            \code{nlm} function, to be passed to it via its \code{iterlim}
            argument.}
  \item{object, x}{an object of class \code{sem} returned by the \code{sem} function.}
  \item{digits}{number of digits for printed output.}
} 

\details{
  The model is set up using RAM ('reticular action model' -- don't ask!) 
  notation -- a simple format
  for specifying general structural equation models by coding the
  `arrows' in the path diagram for the model (see, e.g., McArdle and McDonald, 1984).
  
  The variables in the \eqn{v} vector in the model (typically, the observed and
  unobserved variables, but not error variables) are numbered from 1 to \eqn{m}.  
  the RAM matrix contains one row for each (free or constrained) parameter of the model, and
  five columns, as follows:
    
    \describe{
        \item{1. Number of arrow heads:}{1 (directed arrow) or 2 (covariance).}
        \item{2. Arrow \emph{to}:}{index of the variable at the head of
            a directional arrow, or at one end of a bidirectional arrow.}
        \item{3. Arrow \emph{from}:}{the index of the variable at the tail of 
            a directional arrow, or at the other end of a bidirectonal arrow.}
        \item{4. Parameter number:}{free parameters are numbered from 1 to \eqn{t},
            but do not necessarily appear in consecutive order. Fixed parameters are given
            the number 0. Equality contraints are specified by assigning two 
            parameters the same number.}
        \item{5. Value:}{start value for a free parameter, or value of a fixed parameter. If given
            as NA, the program will compute a start value, by a slight modification of the
            method described by McDonald and Hartmann (1992). \emph{Note:} In some circumstances,
            some start values are selected randomly; this might produce small differences in
            the parameter estimates when the program is rerun.}
        }
        
  \code{sem} fits the model by calling the \code{\link[base]{optim}}
  and \code{\link[base]{nlm}} optimizers
  sequentially to minimize the negative log-likelihood for the model. 
  The \code{L-BFGS-B} method is employed by \code{optim}.
  
  The RAM formulation of the general structural equation model is given by the basic equation
  \deqn{v = Av + u}
  where \eqn{v} and \eqn{u} are vectors of random variables (observed or unobserved), and
  the parameter matrix \eqn{A} contains regression coefficients, symbolized by single-headed arrows
  in a path diagram. Another parameter matrix, 
  \deqn{P = E(uu')}
  contains covariances among the elements of \eqn{u} (assuming that the elements of \eqn{u} have zero
  means). Usually \eqn{v} contains endogenous and exogenous observed and unobserved variables, but not
  error variables (see the examples below).
  
  The \code{start.values} function may be called directly, but is usually called by \code{ram}.
}

\value{
  \code{sem} returns an object of class \code{sem}, with the following elements:
  \item{ram}{RAM matrix, including any rows generated for covariances
    among fixed exogenous variables; column 5 includes computed start values.}
  \item{coef}{estimates of free parameters.}
  \item{criterion}{fitting criterion --- negative log-liklihood divided by \eqn{N - 1}.}
  \item{par}{vector of parameter estimates.}
  \item{cov}{estimated asymptotic covariance matrix of parameter estimates.}
  \item{S}{observed covariance matrix.}
  \item{C}{model-reproduced covariance matrix.}
  \item{A}{RAM \eqn{A} matrix.}
  \item{P}{RAM \eqn{P} matrix.}
  \item{n.fix}{number of fixed exogenous variables.}
  \item{n}{number of observed variables.}
  \item{N}{number of observations.}
  \item{m}{number of variables (observed plus unobserved).}
  \item{t}{number of free parameters.}
  \item{par.posn}{indices of free parameters.}
  \item{convergence}{convergence code returned by \code{nlm} (a code > 2
    indicates a problem).}
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
    In Bollen, K . A. and Long, J. S. (eds.)
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

R.DHP <- R.DHP + t(R.DHP) - diag(diag(R.DHP))   # form symmetric correlation matrix
            
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
                
sem.dhp <- sem(R.DHP, ram.dhp, 329, params.dhp, vars.dhp, fixed.x=5:10)

summary(sem.dhp)
    
##    Model Chisquare =  26.697   Df =  15 Pr(>Chisq) = 0.031302
##    Goodness-of-fit index =  0.98439
##    Adjusted goodness-of-fit index =  0.94275 
##    BIC =  -94.782 
##    
##    Normalized Residuals
##       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    -0.8010 -0.1180  0.0000 -0.0120  0.0398  1.5700 
##    
##    Parameter Estimates
##            Estimate Std Error  z value   Pr(>|z|)                     
##    lam21   1.062673  0.091968 11.55488 0.0000e+00  REdAsp <--- RGenAsp
##    lam42   0.929732  0.071152 13.06682 0.0000e+00  FEdAsp <--- FGenAsp
##    gam11   0.161220  0.038487  4.18892 1.4015e-05 RGenAsp <--- RParAsp
##    gam12   0.249647  0.044581  5.59989 1.0724e-08     RGenAsp <--- RIQ
##    gam13   0.218402  0.043476  5.02350 2.5369e-07    RGenAsp <--- RSES
##    gam14   0.071836  0.050335  1.42714 7.6769e-02    RGenAsp <--- FSES
##    gam23   0.061879  0.051738  1.19602 1.1584e-01    FGenAsp <--- RSES
##    gam24   0.228863  0.044495  5.14361 1.3476e-07    FGenAsp <--- FSES
##    gam25   0.349030  0.044551  7.83443 2.3315e-15     FGenAsp <--- FIQ
##    gam26   0.159529  0.040129  3.97540 3.5131e-05 FGenAsp <--- FParAsp
##    beta12  0.184245  0.096209  1.91505 2.7743e-02 RGenAsp <--- FGenAsp
##    beta21  0.235502  0.119744  1.96671 2.4608e-02 FGenAsp <--- RGenAsp
##    theta1  0.412143  0.052211  7.89385 1.4433e-15 ROccAsp <--> ROccAsp
##    theta2  0.336146  0.053323  6.30399 1.4504e-10   REdAsp <--> REdAsp
##    theta3  0.311197  0.046665  6.66880 1.2895e-11 FOccAsp <--> FOccAsp
##    theta4  0.404601  0.046733  8.65773 0.0000e+00   FEdAsp <--> FEdAsp
##    psi11   0.280987  0.046311  6.06743 6.4985e-10 RGenAsp <--> RGenAsp
##    psi22   0.263832  0.044901  5.87586 2.1033e-09 FGenAsp <--> FGenAsp
##    psi12  -0.022620  0.051650 -0.43795 3.3071e-01 RGenAsp <--> FGenAsp

\dontrun{

# -------------------- Wheaton et al. alienation data ----------------------
    

S.wh <- matrix(c(
   11.834,     0,        0,        0,       0,        0,
    6.947,    9.364,     0,        0,       0,        0,
    6.819,    5.091,   12.532,     0,       0,        0,
    4.783,    5.028,    7.495,    9.986,    0,        0,
   -3.839,   -3.889,   -3.841,   -3.625,   9.610,     0,
  -21.899,  -18.831,  -21.748,  -18.775,  35.522,  450.288)
  , 6, 6)
  
S.wh <- S.wh + t(S.wh) - diag(diag(S.wh))

# This is the model in the SAS manual for PROC CALIS: A Recursive SEM with
# latent endogenous and exogenous variables.
# Curiously, both factor loadings for two of the latent variables are fixed.
 
ram.wh <- matrix(c(             
                1, 1, 7, 0, 1,
                1, 2, 7, 0, .833,
                1, 3, 8, 0, 1,
                1, 4, 8, 0, .833,
                1, 5, 9, 0, 1,
                1, 6, 9, 1, NA,
                1, 7, 9, 2, NA,
                1, 8, 7, 3, NA,
                1, 8, 9, 4, NA,
                2, 1, 1, 5, NA,
                2, 2, 2, 6, NA,
                2, 3, 3, 5, NA,
                2, 4, 4, 6, NA,
                2, 5, 5, 7, NA,
                2, 6, 6, 8, NA,
                2, 1, 3, 9, NA,
                2, 2, 4, 9, NA,
                2, 7, 7, 10, NA,
                2, 8, 8, 11, NA,
                2, 9, 9, 12, NA
                ), ncol=5, byrow=T)

par.names.wh <- c('lamb', 'gam1', 'beta', 'gam2', 
            'the1', 'the2', 'the3', 'the4', 'the5', 'psi1', 'psi2', 'phi')
            
var.names.wh <- c('Anomia67','Powerless67','Anomia71','Powerless71','Education',
    'SEI', 'Alienation67','Alienation71','SES')

sem.wh <- sem(S.wh, ram.wh, 932, par.names.wh, var.names.wh)

summary(sem.wh)

##  Model Chisquare =  13.485   Df =  9 Pr(>Chisq) = 0.14186
##  Goodness-of-fit index =  0.99527
##  Adjusted goodness-of-fit index =  0.98896 
##  BIC =  -64.177
##  
##  Normalized Residuals
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -1.26e+00 -1.31e-01  2.09e-05 -2.87e-02  1.14e-01  8.75e-01 
##  
##  Parameter Estimates
##      Estimate Std Error  z value   Pr(>|z|)                               
##  lamb   5.36882  0.433985  12.3710 0.0000e+00                   SEI <--- SES
##  gam1  -0.62994  0.056128 -11.2233 0.0000e+00          Alienation67 <--- SES
##  beta   0.59311  0.046820  12.6679 0.0000e+00 Alienation71 <--- Alienation67
##  gam2  -0.24086  0.055202  -4.3633 6.4058e-06          Alienation71 <--- SES
##  the1   3.60786  0.200588  17.9864 0.0000e+00         Anomia67 <--> Anomia67
##  the2   3.59494  0.165233  21.7567 0.0000e+00   Powerless67 <--> Powerless67
##  the3   2.99366  0.498975   5.9996 9.8888e-10       Education <--> Education
##  the4 259.57813 18.321380  14.1680 0.0000e+00                   SEI <--> SEI
##  the5   0.90580  0.121710   7.4423 4.9516e-14         Anomia67 <--> Anomia71
##  psi1   5.67049  0.422905  13.4084 0.0000e+00 Alienation67 <--> Alienation67
##  psi2   4.51480  0.334991  13.4774 0.0000e+00 Alienation71 <--> Alienation71
##  phi    6.61633  0.639515  10.3459 0.0000e+00                   SES <--> SES


# The same model, but treating one loading for each latent variable as free.

ram.wh.1 <- matrix(c(           
                1, 1, 7, 0, 1,
                1, 2, 7,13, NA,
                1, 3, 8, 0, 1,
                1, 4, 8,13, NA,
                1, 5, 9, 0, 1,
                1, 6, 9, 1, NA,
                1, 7, 9, 2, NA,
                1, 8, 7, 3, NA,
                1, 8, 9, 4, NA,
                2, 1, 1, 5, NA,
                2, 2, 2, 6, NA,
                2, 3, 3, 5, NA,
                2, 4, 4, 6, NA,
                2, 5, 5, 7, NA,
                2, 6, 6, 8, NA,
                2, 1, 3, 9, NA,
                2, 2, 4, 9, NA,
                2, 7, 7, 10,NA,
                2, 8, 8, 11,NA,
                2, 9, 9, 12,NA
                ), ncol=5, byrow=T)

par.names.wh.1 <- c('lambx', 'gam1', 'beta', 'gam2', 
            'the1', 'the2', 'the3', 'the4', 'the5', 'psi1', 'psi2', 'phi',
            'lamby')

sem.wh.1 <- sem(S.wh, ram.wh.1, 932,  par.names.wh.1, var.names.wh)

summary(sem.wh.1)

##  Model Chisquare =  12.673   Df =  8 Pr(>Chisq) = 0.12360
##  Goodness-of-fit index =  0.99553
##  Adjusted goodness-of-fit index =  0.98828 
##  BIC =  -56.36 
##  
##  Normalized Residuals
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -9.97e-01 -1.40e-01  3.22e-05 -2.87e-02  1.00e-01  7.59e-01 
##  
##  Parameter Estimates
##      Estimate Std Error  z value   Pr(>|z|)                               
##  lambx   5.35309  0.432585  12.3746 0.0000e+00                   SEI <--- SES
##  gam1   -0.62130  0.056143 -11.0665 0.0000e+00          Alienation67 <--- SES
##  beta    0.59427  0.047041  12.6331 0.0000e+00 Alienation71 <--- Alienation67
##  gam2   -0.23581  0.054685  -4.3122 8.0802e-06          Alienation71 <--- SES
##  the1    3.74499  0.249822  14.9907 0.0000e+00         Anomia67 <--> Anomia67
##  the2    3.49380  0.200754  17.4034 0.0000e+00   Powerless67 <--> Powerless67
##  the3    2.97422  0.499650   5.9526 1.3196e-09       Education <--> Education
##  the4  260.12995 18.297861  14.2164 0.0000e+00                   SEI <--> SEI
##  the5    0.90380  0.121819   7.4192 5.8842e-14         Anomia67 <--> Anomia71
##  psi1    5.47364  0.464062  11.7951 0.0000e+00 Alienation67 <--> Alienation67
##  psi2    4.36409  0.362721  12.0315 0.0000e+00 Alienation71 <--> Alienation71
##  phi     6.63578  0.640430  10.3614 0.0000e+00                   SES <--> SES
##  lamby   0.86262  0.033383  25.8400 0.0000e+00  Powerless67 <--- Alienation67 


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
R.thur <- R.thur + t(R.thur) - diag(diag(R.thur))

ram.thur <- matrix(c(
#       heads   to      from    param   start
        1,       1,     10,      1,     NA,     # lam11
        1,       2,     10,      2,     NA,     # lam21
        1,       3,     10,      3,     NA,     # lam31
        1,       4,     11,      4,     NA,     # lam42
        1,       5,     11,      5,     NA,     # lam52
        1,       6,     11,      6,     NA,     # lam62
        1,       7,     12,      7,     NA,     # lam73
        1,       8,     12,      8,     NA,     # lam83
        1,       9,     12,      9,     NA,     # lam93
        1,      10,     13,     10,     NA,     # gam1
        1,      11,     13,     11,     NA,     # gam2
        1,      12,     13,     12,     NA,     # gam3
        2,       1,      1,     13,     NA,     # th1
        2,       2,      2,     14,     NA,     # th2
        2,       3,      3,     15,     NA,     # th3
        2,       4,      4,     16,     NA,     # th4
        2,       5,      5,     17,     NA,     # th5
        2,       6,      6,     18,     NA,     # th6
        2,       7,      7,     19,     NA,     # th7
        2,       8,      8,     20,     NA,     # th8
        2,       9,      9,     21,     NA,     # th9
        2,      10,     10,      0,     1,
        2,      11,     11,      0,     1,
        2,      12,     12,      0,     1,
        2,      13,     13,      0,     1
        ), ncol=5, byrow=T)
        
par.names.thur <- c('lam11', 'lam21', 'lam31', 'lam42', 'lam52', 'lam62', 
                'lam73', 'lam83', 'lam93', 'gam1', 'gam2', 'gam3',
                'th1', 'th2', 'th3', 'th4', 'th5', 'th6', 'th7', 'th8', 'th9')
                
var.names.thur <- c('Sentences','Vocabulary','Sent.Completion','First.Letters',
                '4.Letter.Words','Sufixes','Letter.Series','Pedigrees',
                'Letter.Group', 'F1', 'F2','F3','F4')

        
sem.thur <- sem(R.thur, ram.thur, 213, par.names.thur, var.names.thur)

summary(sem.thur)

##  Model Chisquare =  38.196   Df =  24 Pr(>Chisq) = 0.033101
##  Goodness-of-fit index =  0.95957
##  Adjusted goodness-of-fit index =  0.9242 
##  BIC =  -143.21
##  
##  Normalized Residuals
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -9.75e-01 -4.17e-01  2.51e-05  4.02e-02  9.41e-02  1.63e+00 
##  
##  Parameter Estimates
##      Estimate Std Error z value   Pr(>|z|)                                     
##  lam11  0.51511  0.064966  7.9289 1.1102e-15                    Sentences <--- F1
##  lam21  0.52030  0.065163  7.9845 6.6613e-16                   Vocabulary <--- F1
##  lam31  0.48742  0.062423  7.8083 2.8866e-15              Sent.Completion <--- F1
##  lam42  0.52112  0.063136  8.2540 1.1102e-16                First.Letters <--- F2
##  lam52  0.49707  0.059673  8.3299 0.0000e+00               4.Letter.Words <--- F2
##  lam62  0.43807  0.056479  7.7563 4.3299e-15                      Sufixes <--- F2
##  lam73  0.45243  0.071371  6.3391 1.1552e-10                Letter.Series <--- F3
##  lam83  0.41729  0.061037  6.8367 4.0519e-12                    Pedigrees <--- F3
##  lam93  0.40763  0.064524  6.3175 1.3294e-10                 Letter.Group <--- F3
##  gam1   1.44386  0.264190  5.4652 2.3117e-08                           F1 <--- F4
##  gam2   1.25382  0.216593  5.7888 3.5439e-09                           F2 <--- F4
##  gam3   1.40656  0.279333  5.0354 2.3842e-07                           F3 <--- F4
##  th1    0.18150  0.028400  6.3907 8.2566e-11             Sentences <--> Sentences
##  th2    0.16493  0.027797  5.9334 1.4837e-09           Vocabulary <--> Vocabulary
##  th3    0.26713  0.033468  7.9817 7.7716e-16 Sent.Completion <--> Sent.Completion
##  th4    0.30150  0.050686  5.9484 1.3536e-09     First.Letters <--> First.Letters
##  th5    0.36450  0.052358  6.9617 1.6808e-12   4.Letter.Words <--> 4.Letter.Words
##  th6    0.50641  0.059962  8.4455 0.0000e+00                 Sufixes <--> Sufixes
##  th7    0.39033  0.061599  6.3367 1.1735e-10     Letter.Series <--> Letter.Series
##  th8    0.48137  0.065388  7.3617 9.0816e-14             Pedigrees <--> Pedigrees
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

R.kerch <- R.kerch + t(R.kerch) - diag(diag(R.kerch))
    
ram.kerch<-matrix(c(
#   heads   to      from    param   start
    1,      5,      1,       1,     NA,     # gam51
    1,      5,      2,       2,     NA,     # gam52
    1,      5,      3,       3,     NA,     # gam53
    1,      5,      4,       4,     NA,     # gam54
    1,      6,      1,       5,     NA,     # gam61
    1,      6,      2,       6,     NA,     # gam62
    1,      6,      3,       7,     NA,     # gam63
    1,      6,      4,       8,     NA,     # gam64
    1,      6,      5,       9,     NA,     # beta65
    1,      7,      1,      10,     NA,     # gam71
    1,      7,      2,      11,     NA,     # gam72
    1,      7,      3,      12,     NA,     # gam73
    1,      7,      4,      13,     NA,     # gam74
    1,      7,      5,      14,     NA,     # beta75
    1,      7,      6,      15,     NA,     # beta76
    2,      5,      5,      16,     NA,     # psi5
    2,      6,      6,      17,     NA,     # psi6
    2,      7,      7,      18,     NA),    # psi7
    ncol=5, byrow=T)
    
par.names.kerch <- c('gam51','gam52','gam53','gam54','gam61','gam62',
    'gam63','gam64','beta65','gam71','gam72','gam73','gam74',
    'beta75','beta76', 'psi5','psi6','psi7')
    
var.names.kerch <- c('Intelligence','Siblings','FatherEd','FatherOcc',
    'Grades','EducExp','Occup.Asp')

sem.kerch <- sem(R.kerch, ram.kerch, 737, par.names.kerch, 
    var.names.kerch,fixed.x=1:4)
    
summary(sem.kerch)

##  Model Chisquare =  6.537e-13   Df =  0 Pr(>Chisq) = NA
##  Goodness-of-fit index =  1
##  Adjusted goodness-of-fit index =  NA 
##  BIC =  NA 
##  
##  Normalized Residuals
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  -3.61e-07 -9.77e-14 -7.45e-15 -3.27e-08  1.12e-15  1.27e-13 
##  
##  Parameter Estimates
##          Estimate Std Error  z value   Pr(>|z|)                            
##  gam51   0.525902  0.031182 16.86530 0.0000e+00    Grades <--- Intelligence
##  gam52  -0.029942  0.030149 -0.99314 1.6032e-01        Grades <--- Siblings
##  gam53   0.118966  0.038259  3.10951 9.3699e-04        Grades <--- FatherEd
##  gam54   0.040603  0.037785  1.07456 1.4129e-01       Grades <--- FatherOcc
##  gam61   0.160270  0.032710  4.89979 4.7970e-07   EducExp <--- Intelligence
##  gam62  -0.111779  0.026876 -4.15899 1.5983e-05       EducExp <--- Siblings
##  gam63   0.172719  0.034306  5.03461 2.3941e-07       EducExp <--- FatherEd
##  gam64   0.151852  0.033688  4.50758 3.2785e-06      EducExp <--- FatherOcc
##  beta65  0.405150  0.032838 12.33799 0.0000e+00         EducExp <--- Grades
##  gam71  -0.039405  0.034500 -1.14215 1.2670e-01 Occup.Asp <--- Intelligence
##  gam72  -0.018825  0.028222 -0.66700 2.5238e-01     Occup.Asp <--- Siblings
##  gam73  -0.041333  0.036216 -1.14126 1.2688e-01     Occup.Asp <--- FatherEd
##  gam74   0.099577  0.035446  2.80924 2.4829e-03    Occup.Asp <--- FatherOcc
##  beta75  0.157912  0.037443  4.21738 1.2358e-05       Occup.Asp <--- Grades
##  beta76  0.549593  0.038260 14.36486 0.0000e+00      Occup.Asp <--- EducExp
##  psi5    0.650995  0.033946 19.17743 0.0000e+00          Grades <--> Grades
##  psi6    0.516652  0.026943 19.17590 0.0000e+00        EducExp <--> EducExp
##  psi7    0.556617  0.029026 19.17644 0.0000e+00    Occup.Asp <--> Occup.Asp 

    }
}

\keyword{models}
