# last modified 14 April 2001 by J. Fox

residuals.sem <- function(object, ...){
    object$S - object$C
    }    

standardized.residuals <- function(object, ...){
    UseMethod("standardized.residuals")
    }

standardized.residuals.sem <- function(object, ...){
    res <- residuals(object)
    s <- diag(object$S)
    res/sqrt(outer(s, s))
    }

normalized.residuals <- function(object, ...){
    UseMethod("normalized.residuals")
    }
    
normalized.residuals.sem <- function(object, ...){
    res <- residuals(object)
    N <- object$N
    C <- object$C
    c <- diag(C)
    res/sqrt((outer(c,c) + C^2)/N)
    }
