# last modified 29 Jan 2002 by J. Fox

print.sem <- function(x, ...) {
    n <- x$n
    t <- x$t
    n.fix <- x$n.fix
    df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
    cat("\n Model Chisquare = ", x$criterion * (x$N - 1), 
        "  Df = ", df, "\n\n")
    print(x$coeff)
    cat("\n Iterations = ", x$iterations, "\n")
    invisible(x)
    }
