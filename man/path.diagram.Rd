\name{path.diagram}

\alias{path.diagram}
\alias{path.diagram.sem}

\title{Draw Path Diagram}
\description{
  \code{path.diagram} creates a description of the path diagram
  for a structural-equation-model object to be processed by the
  graph-drawing program \emph{dot}; see Koutsofios and North (1996)
  and \url{http://www.research.att.com/sw/tools/graphviz/}.
}

\usage{
path.diagram(model, ...)

path.diagram.sem(model, out.file, min.rank=NULL, max.rank=NULL, 
    same.rank=NULL, variables=model$var.names, parameters=rownames(object$ram), 
    ignore.double=T, edge.labels=c("names", "values"), size=c(8, 8), 
    node.font=c("Helvetica", 14), edge.font=c("Helvetica", 10), 
    rank.direction=c("LR", "TB"))
}

\arguments{
  \item{model}{a structural-equation-model object produced by \code{sem}.}
  \item{\dots}{arguments to be passed down to \code{path.diagram.sem}.}
  \item{out.file}{a file to which to write the \emph{dot} description of the
    path diagram; if not specified, the description is written to standard
    output (normally the R Console).}
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
    \code{"values"} to label arrows with parameter estimates.}
  \item{size}{the size of the graph, in inches.}
  \item{node.font}{font name and point-size for printing variable names.}
  \item{edge.font}{font name and point-size for printing arrow names or values.}
  \item{rank.direction}{draw graph left-to-right, \code{"LR"}, the default, 
    or top-to-bottom, \code{"TB"}.}
}

}
\value{
  \code{NULL}: \code{path.diagram} is used for its side-effect, producing a
  graph description of the model.}
}

\references{
    Koutsofios, E., and North, S. C. (1996)
    Drawing graphs with \emph{dot}.
    \url{http://www.research.att.com/sw/tools/graphviz/dotguide.pdf}.
    
}
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{

# The Duncan, Haller, and Portes Peer-Influences Model

R.DHP <- matrix(c(     
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

rownames(R.DHP) <- colnames(R.DHP) <- c('ROccAsp', 'REdAsp', 'FOccAsp', 
    'FEdAsp', 'RParAsp', 'RIQ', 'RSES', 'FSES', 'FIQ', 'FParAsp')

sem.dhp <- sem(model.dhp, R.DHP, 329,
    fixed.x=c('RParAsp', 'RIQ', 'RSES', 'FSES', 'FIQ', 'FParAsp'))
    
path.diagram(sem.dhp, min.rank='RIQ, RSES, RParAsp, FParAsp, FSES, FIQ', 
    max.rank='ROccAsp, REdAsp, FEdAsp, FOccAsp')

##   digraph "sem.dhp" {
##   rankdir=LR;
##   size="8,8";
##   node [fontname="Helvetica" fontsize=14 shape=box];
##   edge [fontname="Helvetica" fontsize=10];
##   center=1;
##   {rank=min "RIQ" "RSES" "RParAsp" "FParAsp" "FSES" "FIQ"}
##   {rank=max "ROccAsp" "REdAsp" "FEdAsp" "FOccAsp"}
##   "RGenAsp" [shape=ellipse]
##   "FGenAsp" [shape=ellipse]
##   "RParAsp" -> "RGenAsp" [label="gam11"];
##   "RIQ" -> "RGenAsp" [label="gam12"];
##   "RSES" -> "RGenAsp" [label="gam13"];
##   "FSES" -> "RGenAsp" [label="gam14"];
##   "RSES" -> "FGenAsp" [label="gam23"];
##   "FSES" -> "FGenAsp" [label="gam24"];
##   "FIQ" -> "FGenAsp" [label="gam25"];
##   "FParAsp" -> "FGenAsp" [label="gam26"];
##   "FGenAsp" -> "RGenAsp" [label="beta12"];
##   "RGenAsp" -> "FGenAsp" [label="beta21"];
##   "RGenAsp" -> "ROccAsp" [label=""];
##   "RGenAsp" -> "REdAsp" [label="lam21"];
##   "FGenAsp" -> "FOccAsp" [label=""];
##   "FGenAsp" -> "FEdAsp" [label="lam42"];
##   }

}
\keyword{dplot}
\keyword{models}
