# last modified 14 April 2001 by J. Fox

print.sem <- function(x) {
    n <- x$n
    t <- x$t
    cat("\n Model Chisquare = ", x$criterion * (x$N - 1), 
        "  Df = ", n*(n + 1)/2 - t, "\n\n")
    print(x$coeff)
    invisible(x)
    }
