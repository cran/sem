# last modified 1 July 2007 by J. Fox
                                                               
summary.sem <- function(object, digits=5, conf.level=.90, ...) {
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
    dfNull <- n*(n - 1)/2
    invC <- solve(C)
    CSC <- invC %*% (S - C)
    CSC <- CSC %*% CSC
    CS <- invC %*% S
    CS <- CS %*% CS
    chisq <- object$criterion * (N - (!object$raw))
    chisqNull <- object$chisqNull
    GFI <- if (!object$raw) 1 - sum(diag(CSC))/sum(diag(CS)) else NA
    if ((!object$raw) && df > 0){
        AGFI <- 1 - (n*(n + 1)/(2*df))*(1 - GFI)
        NFI <- (chisqNull - chisq)/chisqNull
        NNFI <- (chisqNull/dfNull - chisq/df)/(chisqNull/dfNull - 1)
        L1 <- max(chisq - df, 0)
        L0 <- max(L1, chisqNull - dfNull)
        CFI <- 1 - L1/L0
        RMSEA <- sqrt(max(object$criterion/df - 1/(N - (!object$raw)), 0))
        tail <- (1 - conf.level)/2 
        max <- N
        while (max > 1){
            res <- optimize(function(lam) (tail - pchisq(chisq, df, ncp=lam))^2, interval=c(0, max))
            if (sqrt(res$objective) < tail/100) break
            max <- max/2
            }
        lam.U <- if (max <= 1) NA else res$minimum
        max <- max(max, 1)
        while (max > 1){
            res <- optimize(function(lam) (1 - tail - pchisq(chisq, df, ncp=lam))^2, interval=c(0, max))
            if (sqrt(res$objective) < tail/100) break
            max <- max/2
            }
        lam.L <- if (max <= 1) NA else res$minimum
        RMSEA.U <- sqrt(lam.U/((N - (!object$raw))*df))
        RMSEA.L <- sqrt(lam.L/((N - (!object$raw))*df))
        }
    else RMSEA.U <- RMSEA.L <- RMSEA <- NFI <- NNFI <- CFI <- AGFI <- NA
    RMSEA <- c(RMSEA, RMSEA.L, RMSEA.U, conf.level)
    if (!is.null(object$coeff)){
        var.names <- rownames(object$A)
        ram <- object$ram[object$par.posn, , drop=FALSE]
        par.code <- paste(var.names[ram[,2]], c('<---', '<-->')[ram[,1]],
                        var.names[ram[,3]])
        coeff <- data.frame(object$coeff, se, z, 2*(1 - pnorm(abs(z))), par.code)
        names(coeff) <- c("Estimate", "Std Error", "z value", "Pr(>|z|)", " ")
        row.names(coeff) <- names(object$coeff)
        }
    else coeff <- NULL
    BIC <- chisq - df * log(N)
    SRMR <- sqrt(sum(standardized.residuals(object)^2 * 
        upper.tri(diag(n), diag=TRUE))/(n*(n + 1)/2))
    ans <- list(chisq=chisq, df=df, chisqNull=chisqNull, dfNull=dfNull,
        GFI=GFI, AGFI=AGFI, RMSEA=RMSEA, NFI=NFI, NNFI=NNFI, CFI=CFI, BIC=BIC, SRMR=SRMR, 
        norm.res=norm.res, coeff=coeff, digits=digits, 
        iterations=object$iterations, aliased=object$aliased, raw=object$raw)
    class(ans) <- "summary.sem"
    ans
    }
    
print.summary.sem <- function(x, ...){
    old.digits <- options(digits=x$digits)
    on.exit(options(old.digits))
    if (x$raw) cat("\nModel fit to raw moment matrix.\n")
    cat("\n Model Chisquare = ", x$chisq, "  Df = ", x$df, 
        "Pr(>Chisq) =", if (x$df > 0) 1 - pchisq(x$chisq, x$df)
            else NA)
    if (!x$raw) {
       cat("\n Chisquare (null model) = ", x$chisqNull,  "  Df = ", x$dfNull)
       cat("\n Goodness-of-fit index = ", x$GFI)
       }
    if (x$df > 0 && !x$raw){
        cat("\n Adjusted goodness-of-fit index = ", x$AGFI)
        cat("\n RMSEA index =  ", x$RMSEA[1],
            "   ", 100*x$RMSEA[4], "% CI: (", x$RMSEA[2], ", ", x$RMSEA[3],")", sep="")
        cat("\n Bentler-Bonnett NFI = ", x$NFI)
        cat("\n Tucker-Lewis NNFI = ", x$NNFI)
        cat("\n Bentler CFI = ", x$CFI)
        cat("\n SRMR = ", x$SRMR)
        }
    cat("\n BIC = ", x$BIC, "\n")
    cat("\n Normalized Residuals\n")
    print(summary(as.vector(x$norm.res)))
    if (!is.null(x$coeff)){
        cat("\n Parameter Estimates\n")
        print(x$coeff, right=FALSE)
        cat("\n Iterations = ", x$iterations, "\n")
        if (!is.null(x$aliased)) cat("\n Aliased parameters:", x$aliased, "\n")
        }
    invisible(x)
    }
    
deviance.sem <- function(object, ...) object$criterion * (object$N - (!object$raw))

df.residual.sem <- function(object, ...) {
    n.fix <- object$n.fix
    n <- object$n
    t <- object$t
    n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
    }
