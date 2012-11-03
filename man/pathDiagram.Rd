\name{pathDiagram}

\alias{pathDiagram}
\alias{pathDiagram.sem}
\alias{pathDiagram.semmod}

\title{Draw Path Diagram}
\description{
  \code{pathDiagram} creates a description of the path diagram
  for a structural-equation-model or SEM-specification object to be processed by the
  graph-drawing program \emph{dot}, which can be called
  automatically; see Koutsofios and North (2002)
  and \url{http://www.graphviz.org/}. To obtain graphics output
  directly, the \emph{dot} program must be on the system search path. The \code{semmod}
  method sets up a call to the \code{sem} method.
}

\usage{
\method{pathDiagram}{sem}(model, file, min.rank=NULL, max.rank=NULL,
	same.rank=NULL, variables=model$var.names, parameters,
	ignore.double=TRUE, edge.labels=c("names", "values", "both"), 
	size=c(8, 8), node.font=c("Helvetica", 14),
	edge.font=c("Helvetica", 10), rank.direction=c("LR", "TB"), 
	digits=2, standardize=FALSE, output.type=c("graphics", "dot"),
	graphics.fmt="pdf", dot.options=NULL, ...)

\method{pathDiagram}{semmod}(model, obs.variables, ...)
}

\arguments{
  \item{model}{a structural-equation-model or SEM-specification object produced by \code{sem}, 
    or, respectively, \code{specifyEquations}, \code{specifyModel}, or \code{cfa}.}
  \item{\dots}{arguments passed down, e.g., from the \code{semmod} method to the \code{sem} method.}
  \item{file}{a file name, given \emph{without} an extension, to which to write the \emph{dot} description of the
    path diagram and, by default, to use that to draw a graph of the model; 
    if not specified, the \emph{dot} description is written to standard output (normally the R Console).
    Note that by default (see the \code{output.type} argument below), \emph{two} files are written:
    a ".dot" file and a graphics file of type specified by the \code{graphics.fmt} argument (below);
    \code{file} may include a path specification.}
  \item{min.rank}{a character string listing names of variables to be assigned
    minimum rank (order) in the graph; the names should be separated by commas.}
  \item{max.rank}{a character string listing names of variables to be assigned
    maximum rank in the graph; the names should be separated by commas.}
  \item{same.rank}{a character string or vector of character strings of variables
    to be assigned equivalent rank in the graph; names in each string should be
    separated by commas.}
  \item{variables}{variable names; defaults to the variable names in \code{model}.
    If specified, the variable names should be in the same order as in \code{model}.}
  \item{parameters}{parameter names; defaults to the parameter names in 
    \code{model}. If specified, the parameter names should be in the same order 
    as in \code{model}.}
  \item{ignore.double}{if \code{TRUE}, the default, double-headed arrows, representing
    variances and covariances, are not graphed.}
  \item{edge.labels}{\code{"names"} to label arrows with parameter names; 
    \code{"values"} to label arrows with parameter estimates, or \code{"both"}.}
  \item{size}{the size of the graph, in inches.}
  \item{node.font}{font name and point-size for printing variable names.}
  \item{edge.font}{font name and point-size for printing arrow names or values.}
  \item{rank.direction}{draw graph left-to-right, \code{"LR"}, the default, 
    or top-to-bottom, \code{"TB"}.}
  \item{digits}{number of digits after the decimal point (default, 2) to which 
    to round parameter estimates.}
  \item{standardize}{if \code{TRUE}, display standardized coefficients; default is \code{FALSE}.}
  \item{output.type}{if \code{"graphics"}, the default, \emph{both} a ".dot" file and a
  	graphics file will be created.}
  \item{graphics.fmt}{a graphics format recognized by the \emph{dot} program; the default is \code{"pdf"};
  	\code{graphics.fmt} is also used for the extension of the graphics file that is created.}
  \item{dot.options}{options to be passed to the \emph{dot} program, given as a character string.}
  \item{obs.variables}{a character vector with the names of the observed variables in the model.}
}

\value{
  \code{NULL}: \code{pathDiagram} is used for its side-effect, producing a
  graph description of the model.
}

\references{
    Koutsofios, E., and North, S. C. (2002)
    Drawing graphs with \emph{dot}.
    \url{http://www.graphviz.org/Documentation/dotguide.pdf}.    
}

\author{John Fox \email{jfox@mcmaster.ca}, Adam Kramer, and Michael Friendly}

\seealso{\code{\link{sem}}, \code{\link{specifyEquations}}, \code{\link{specifyModel}}, \code{\link{cfa}}}

\examples{
	\dontrun{
# The Duncan, Haller, and Portes Peer-Influences Model

R.DHP <- readMoments(diag=FALSE, names=c('ROccAsp', 'REdAsp', 'FOccAsp', 
                'FEdAsp', 'RParAsp', 'RIQ', 'RSES', 'FSES', 'FIQ', 'FParAsp'))
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
    RGenAsp  -> ROccAsp,  NA,       1
    RGenAsp  -> REdAsp,  lam21,  NA
    FGenAsp  -> FOccAsp,  NA,       1
    FGenAsp  -> FEdAsp,  lam42,  NA
    RGenAsp <-> RGenAsp, ps11,   NA
    FGenAsp <-> FGenAsp, ps22,   NA
    RGenAsp <-> FGenAsp, ps12,   NA
    ROccAsp <-> ROccAsp, theta1, NA
    REdAsp  <-> REdAsp,  theta2, NA
    FOccAsp <-> FOccAsp, theta3, NA
    FEdAsp  <-> FEdAsp,  theta4, NA

sem.dhp <- sem(model.dhp, R.DHP, 329,
    fixed.x=c('RParAsp', 'RIQ', 'RSES', 'FSES', 'FIQ', 'FParAsp'))
    
pathDiagram(sem.dhp, min.rank='RIQ, RSES, RParAsp, FParAsp, FSES, FIQ', 
    max.rank='ROccAsp, REdAsp, FEdAsp, FOccAsp')
    
      # equivalent:    
pathDiagram(model.dhp, 
    obs.variables=c("RParAsp", "RIQ", "RSES", "FSES", "FIQ", 
      "FParAsp", "ROccAsp", "REdAsp", "FOccAsp", "FEdAsp"),
    min.rank='RIQ, RSES, RParAsp, FParAsp, FSES, FIQ', 
    max.rank='ROccAsp, REdAsp, FEdAsp, FOccAsp',
    same.rank='RGenAsp, FGenAsp') # note crowding of labels
    }
}

\keyword{dplot}
\keyword{models}
