# last modified 10 June 05 by J. Fox

raw.moments <- function(object, ...) UseMethod("raw.moments")

raw.moments.formula <- function (formula, data, subset, na.action, 
        contrasts = NULL, ...) {
    if (missing(na.action))
        na.action <- options()$na.action
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, sys.frame(sys.parent()))))
        m$data <- as.data.frame(data)
    m$instruments <- m$contrasts <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, sys.frame(sys.parent()))
    response <- attr(attr(mf, "terms"), "response")
    if (response) stop("formula cannot have a response")
    na.act <- attr(mf, "na.action")
    raw.moments(model.matrix(formula, data = mf, contrasts))
    }
    
raw.moments.default <- function(object, ...){
    object <- as.matrix(object)
    N <- nrow(object)
    result <- crossprod(object, object)/N
    attr(result, "N") <- N
    class(result) <- "rawmoments"
    result
    }

print.rawmoments <- function(x, ...){
    xx <- unclass(x)
    attr(xx, "N") <- NULL
    cat("\nRaw Moments\n")
    print(xx)
    cat("\nN = ", attr(x, "N"), "\n")
    invisible(x)
    }
    

