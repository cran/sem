\name{pathDiagram}

\alias{pathDiagram}
\alias{pathDiagram.sem}
\alias{pathDiagram.semmod}
\alias{math}

\title{Draw Path Diagram}
\description{
  \code{pathDiagram} creates a description of the path diagram
  for a structural-equation-model or SEM-specification object to be processed by the
  graph-drawing program \emph{dot}.
}

\usage{
pathDiagram(model, ...)

\method{pathDiagram}{sem}(model, file = "pathDiagram", 
    style = c("ram", "traditional"),
    output.type = c("html", "graphics", "dot"), graphics.fmt = "pdf", 
    dot.options = NULL,
    size = c(8, 8), node.font = c("Helvetica", 14),
    edge.font = c("Helvetica", 10), digits = 2, 
    rank.direction = c("LR", "TB"),
    min.rank = NULL, max.rank = NULL, same.rank = NULL,
    variables = model$var.names, var.labels, parameters, par.labels,
    ignore.double = TRUE, ignore.self = FALSE, error.nodes = TRUE,
    edge.labels = c("names", "values", "both"),  
    edge.colors = c("black", "black"),
    edge.weight = c("fixed", "proportional"),
    node.colors = c("transparent", "transparent", "transparent"),
    standardize = FALSE, ...)

\method{pathDiagram}{semmod}(model, obs.variables, ...)

math(text, html.only=FALSE, hat=FALSE)
}

\arguments{
  \item{model}{a structural-equation-model or SEM-specification object produced by \code{sem}, 
    or, respectively, \code{specifyEquations}, \code{specifyModel}, or \code{cfa}.}
  \item{\dots}{arguments passed down, e.g., from the \code{semmod} method to the \code{sem} method.}
  \item{file}{a file name, by default \code{"pathDiagram"}, 
    given \emph{without} an extension, to which to write the \emph{dot} description of the
    path diagram if \code{output.type} \code{"graphics"} or \code{"dot"} is selected, and for the
    graphics output file (with appropriate extension) if \code{"graphics"} output is selected, in which case
    a "\code{.dot}" file and a graphics file of type specified by the \code{graphics.fmt} argument (below);
    \code{file} may include a path specification.}
  \item{style}{\code{"ram"} (the default) for a RAM path diagram including self-directed double-headed arrows
    representing variances, including error variances; or \code{"traditional"} for a path diagram 
    including nodes representing error variables.}
  \item{output.type}{if \code{"html"} (the default), the path diagram will open in the user"s default 
    web browser; if \code{"dot"}, a file containing \emph{dot} commands will be written; if \code{"graphics"},
    both \code{.dot} and graphics files will be written.}
  \item{graphics.fmt}{a graphics format recognized by the \emph{dot} program; the default is \code{"pdf"};
  	\code{graphics.fmt} is also used for the extension of the graphics file that is created.}
  \item{dot.options}{options to be passed to the \emph{dot} program, given as a character string.}
  \item{size}{the size of the graph, in inches.}
  \item{node.font}{font name and point-size for printing variable names.}
  \item{edge.font}{font name and point-size for printing arrow names or values.}
  \item{digits}{number of digits after the decimal point (default, 2) to which 
    to round parameter estimates.}
  \item{rank.direction}{draw graph left-to-right, \code{"LR"}, the default, 
    or top-to-bottom, \code{"TB"}.}
  \item{min.rank}{a character string listing names of variables to be assigned
    minimum rank (order) in the graph; the names should be separated by commas.}
  \item{max.rank}{a character string listing names of variables to be assigned
    maximum rank in the graph; the names should be separated by commas.}
  \item{same.rank}{a character string or vector of character strings of variables
    to be assigned equivalent rank in the graph; names in each string should be
    separated by commas.}
  \item{variables}{variable names; defaults to the variable names in \code{model}.
    If specified, the variable names should be in the same order as in \code{model}.}
  \item{var.labels}{a character vector with labels to be used 
    in lieu of (some of) the variables names, for greater flexibility
    in labelling nodes in the graph --- e.g., the labels can be created with the \code{math}
    function. The elements of the vector must have names corresponding to variables in the model.}
  \item{parameters}{parameter names; defaults to the parameter names in 
    \code{model}. If specified, the parameter names should be in the same order 
    as in \code{model}.}
  \item{par.labels}{a character vector with labels to be used 
    in lieu of (some of) the parameter names, for greater flexibility
    in labelling edges in the graph --- e.g., the labels can be created with the \code{math}
    function. The elements of the vector must have names corresponding to parameters in the model.}
  \item{ignore.double}{if \code{TRUE}, the default, double-headed arrows, representing
    variances and covariances, are not graphed.}
  \item{ignore.self}{if \code{TRUE} (the default is \code{FALSE}), and \code{ignore.double=FALSE},
    self-directed double-headed arrows representing error variances are suppressed; note that
    if \code{ignore.double=TRUE}, \emph{all} double-headed arrows, including self-directed arrows,
    are suppressed.}
  \item{error.nodes}{if \code{TRUE} (the default) and \code{style="traditional"}, show the nodes
    representing error variables.}
  \item{edge.labels}{\code{"names"} to label arrows with parameter names; 
    \code{"values"} to label arrows with parameter estimates, or \code{"both"}.}
  \item{edge.colors}{two-element character vector giving colors of positive
    and negative arrows respectively; the default is \code{c("black", "black")}.}
  \item{edge.weight}{if \code{"proportional"} (the default is \code{"fixed"}), the thickness of
    edges is proportional to the absolute size of the corresponding parameter estimate;
    this is generally sensible only if \code{standardize=TRUE}.}
  \item{node.colors}{a two- or three-element character vector giving colors of nodes representing
    exogenous, endogenous, and error variables (for traditional path diagrams) consecutively; 
    the default is \code{"transparent"} for all three; if a two colors are given, error variables
    are colored as exogenous (the first color.}
  \item{standardize}{if \code{TRUE}, display standardized coefficients; default is \code{FALSE}.}

  \item{obs.variables}{a character vector with the names of the observed variables in the model.}
  
  \item{text}{a character string or vector of character strings to be translated into node or edge label symbols.
    If a vector of character strings is supplied, then the elements of the vector should be named with the
    corresponding variable (node) or parameter (edge) name.}
  \item{html.only}{If \code{TRUE} (the default is \code{FALSE}), the
    character strings in \code{text} are to be treated as an HTML character codes, in
    which case the prefix \code{"#"} and suffix \code{";"} are appended to each. Otherwise, \code{text} should
    only contain the names of lowercase or uppercase Greek letters, such as \code{"alpha"} or \code{"Alpha"}.
    The full set of Greek letters recognized is given in the file \code{Greek.txt} in the package's
    \code{etc} subdirectory -- or type \code{sem:::Greek} at the R command prompt.
    In either case, the symbols may be followed by numeric subscripts in curly braces consisting of numerals 
    (e.g., \code{"beta_{12}"}),
    and/or numeric superscripts (e.g., \code{"sigma^{2}"},
    \code{"sigma_{1}^{2}"}). Depending upon your OS, subscripts and superscripts may only work properly with HTML output
    from \code{pathDiagram}, not with graphics output produced by \emph{dot}.}
  \item{hat}{If \code{TRUE} (the default is \code{FALSE}), a hat (circumflex) is placed over the
    symbols in \code{text}; this feature doesn't produce a visually appealing result.}
}

\details{
\code{pathDiagram} creates a description of the path diagram
  for a structural-equation-model or SEM-specification object to be processed by the
  graph-drawing program \emph{dot}, which can be called
  automatically; see Koutsofios and North (2002)
  and \url{http://www.graphviz.org/}. To obtain graphics output
  directly, the \emph{dot} program must be on the system search path. 
  
  Alternatively, \emph{HTML} output can be created in a web browser without an independent installation
  of \code{dot} using facilities in the \pkg{DiagrammeR} package.
  
  The \code{math} function can be used to create node (variable) and edge (arrow) labels with
  symbols such as Greek letters, subscripts, and superscripts.
  
  The \code{semmod} method of \code{pathDiagram} sets up a call to the \code{sem} method. 
  
  The various
  arguments to \code{pathDiagram} can be used to customize the diagram, but if there are too many constraints
  on node placement, \emph{dot} may fail to produce a graph or may produce a distorted graph. 
  \code{pathDiagram} can create both RAM-style diagrams, in which variances are represented as self-directed
  arrows, and traditional path diagrams, in which error variables appear explicitly as nodes. As is conventional,
  latent variables (including error variables) are represented as ellipses and observed variables as 
  rectangles; double-headed arrows represent covariances (and in RAM diagrams, variances) and single-headed
  arrows represent structural coefficients.
}

\value{
  \code{pathDiagram} invisibly returns a character vector containing \emph{dot} commands.
  \code{math} returns a character vector containing suitable \code{HTML} markup.
}

\references{
    Koutsofios, E., and North, S. C. (2002)
    Drawing graphs with \emph{dot}.
    \url{http://www.graphviz.org/Documentation/dotguide.pdf}.    
}

\author{John Fox \email{jfox@mcmaster.ca}, Adam Kramer, and Michael Friendly}

\seealso{\code{\link{sem}}, \code{\link{specifyEquations}}, \code{\link{specifyModel}}, \code{\link{cfa}}}

\examples{
	if (interactive()) {
# The Duncan, Haller, and Portes Peer-Influences Model

R.DHP <- readMoments(diag=FALSE, names=c("ROccAsp", "REdAsp", "FOccAsp", 
                "FEdAsp", "RParAsp", "RIQ", "RSES", "FSES", "FIQ", "FParAsp"),
                text="
    .6247                                                              
    .3269  .3669                                                        
    .4216  .3275  .6404                                      
    .2137  .2742  .1124  .0839                                
    .4105  .4043  .2903  .2598  .1839                          
    .3240  .4047  .3054  .2786  .0489  .2220                    
    .2930  .2407  .4105  .3607  .0186  .1861  .2707              
    .2995  .2863  .5191  .5007  .0782  .3355  .2302  .2950        
    .0760  .0702  .2784  .1988  .1147  .1021  .0931 -.0438  .2087  
")

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
")

sem.dhp <- sem(model.dhp, R.DHP, 329,
    fixed.x=c("RParAsp", "RIQ", "RSES", "FSES", "FIQ", "FParAsp"))
    
pathDiagram(sem.dhp, min.rank="RIQ, RSES, RParAsp, FParAsp, FSES, FIQ", 
                     max.rank="ROccAsp, REdAsp, FEdAsp, FOccAsp",
                     same.rank="RGenAsp, FGenAsp",
                     edge.labels="values")
    
pathDiagram(model.dhp,
    obs.variables=c("RParAsp", "RIQ", "RSES", "FSES", "FIQ", 
        "FParAsp", "ROccAsp", "REdAsp", "FOccAsp", "FEdAsp"),
    style="traditional", 
    node.colors=c("pink", "lightblue", "lightgreen"),
    min.rank="RIQ, RSES, RParAsp, FParAsp, FSES, FIQ",
    max.rank="ROccAsp, REdAsp, FEdAsp, FOccAsp",
    same.rank="RGenAsp, FGenAsp",
    var.labels=c(RParAsp="Respondent Parental Aspiration", 
        RIQ="Respondent IQ",
        RSES="Respondent SES",
        FSES="Friend SES",
        FIQ="Friend IQ",
        FParAsp="Friend Parental Aspiration",
        ROccAsp="Respondent Occupational Aspiration",
        REdAsp="Respondent Educational Aspiration",
        RGenAsp="Respondent General Aspiration",
        FOccAsp="Friend Occupational Aspiration",
        FEdAsp="Friend Educational Aspiration",
        FGenAsp="Friend General Aspiration",
        math(c(RGenAsp.error="xi_{1}",
        FGenAsp.error="xi_{2}",
        ROccAsp.error="epsilon_{1}",
        REdAsp.error="epsilon_{2}",
        FOccAsp.error="epsilon_{3}",
        FEdAsp.error="epsilon_{4}"))),
    par.labels=math(c(gam11="gamma_{11}",
        gam12="gamma_{12}",
        gam13="gamma_{13}",
        gam14="gamma_{14}",
        gam23="gamma_{23}",
        gam24="gamma_{24}",
        gam25="gamma_{25}",
        gam26="gamma_{26}",
        beta12="beta_{12}",
        beta21="beta_{21}",
        lam21="lambda_{21}",
        lam42="lambda_{42}",
        ps11="psi_{11}",
        ps22="psi_{22}",
        ps12="psi_{12}",
        theta1="theta_{1}",
        theta2="theta_{2}",
        theta3="theta_{3}",
        theta4="theta_{4}")))
        
    # the following example contributed by Michael Friendly:
    
union <- readMoments(diag=TRUE,
    names=c('y1', 'y2', 'y3', 'x1', 'x2'), text="
14.610
-5.250  11.017
-8.057  11.087   31.971
-0.482   0.677    1.559   1.021
-18.857  17.861   28.250   7.139  215.662
")

union.mod <- specifyEquations(covs=c("x1, x2"), text="
y1 = gam12*x2
y2 = beta21*y1 + gam22*x2
y3 = beta31*y1 + beta32*y2 + gam31*x1
")

union.sem <- sem(union.mod, union, N=173)

dot <- pathDiagram(union.sem, style="traditional",  
    ignore.double=FALSE, error.nodes=FALSE,
    edge.labels="values", 
    min.rank=c("Years", "Age"), 
    max.rank=c("Sentiment", "Sentiment.error"),
    same.rank=c("Deference, Deference.error", "Activism, Activism.error"),
    variables=c("Deference", "Activism", "Sentiment", "Years", "Age"),
    edge.colors=c("black", "red"),
    node.colors = c("pink", "lightblue"))

cat(paste(dot, collapse="\n")) # dot commands

    }
}

\keyword{dplot}
\keyword{models}
