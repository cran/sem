# last modified 2 June 2005 by J. Fox

print.sem <- function(x, ...) {
    n <- x$n
    t <- x$t
    n.fix <- x$n.fix
    df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
    cat("\n Model Chisquare = ", x$criterion * (x$N - (!x$raw)), 
        "  Df = ", df, "\n\n")
    print(x$coeff)
    cat("\n Iterations = ", x$iterations, "\n")
    if (!is.null(x$aliased)) cat("\n Aliased parameters:", x$aliased, "\n")
    invisible(x)
    }
