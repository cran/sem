# last modified 25 July 2001 by J. Fox

summary.sem <- function(object, digits=5) {
    norm.res <- normalized.residuals(object)
    se <- sqrt(diag(object$cov))
    z <- object$coeff/se
    n.fix <- object$n.fix
    n <- object$n
    t <- object$t
    S <- object$S
    C <- object$C
    N <- object$N
    df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
    invC <- solve(C)
    CSC <- invC %*% (S - C)
    CSC <- CSC %*% CSC
    CS <- invC %*% S
    CS <- CS %*% CS
    GFI <- 1 - sum(diag(CSC))/sum(diag(CS))
    AGFI <- if (df > 0) 1 - (n*(n + 1)/(2*df))*(1 - GFI)
     else NA
    var.names <- rownames(object$A)
    ram <- object$ram[object$par.posn,]
    par.code <- paste(var.names[ram[,2]], c('<---', '<-->')[ram[,1]],
                    var.names[ram[,3]])
    coeff <- data.frame(object$coeff, se, z, 1 - pnorm(abs(z)), par.code)
    names(coeff) <- c("Estimate", "Std Error", "z value", "Pr(>|z|)", " ")
    row.names(coeff) <- names(object$coeff)
    chisq <- object$criterion * (N - 1)
    BIC <- if (df > 0) chisq - df * log(N*n) else NA
    ans <- list(chisq=chisq, df=df, GFI=GFI, AGFI=AGFI, BIC=BIC, 
        norm.res=norm.res, coeff=coeff, digits=digits)
    class(ans) <- "summary.sem"
    ans
    }
    
print.summary.sem <- function(x){
    old.digits <- options(digits=x$digits)
    on.exit(options(old.digits))
    cat("\n Model Chisquare = ", x$chisq, "  Df = ", x$df, 
        "Pr(>Chisq) =", if (x$df > 0) 1 - pchisq(x$chisq, x$df)
            else NA)
    cat("\n Goodness-of-fit index = ", x$GFI)
    cat("\n Adjusted goodness-of-fit index = ", x$AGFI)
    cat("\n BIC = ", x$BIC, "\n")
    cat("\n Normalized Residuals\n")
    print(summary(as.vector(x$norm.res)))
    cat("\n Parameter Estimates\n")
    print(x$coeff)
    invisible(x)
    }
