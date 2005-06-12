# bootstrapped standard errors and confidence intervals for sem

# last modified 8 June 2005 by J. Fox

boot.sem <- function(data, model, R=100, cov=cov, ...){
    refit <- function(){
        indices <- sample(N, N, replace=TRUE)
        S <- cov(data[indices,])
        refitted.model <- sem(ram, S, N, coef.names, ...)
        refitted.model$coeff
        }
    if (!require("boot")) stop("package boot not available")
    # the following 2 lines borrowed from boot in package boot
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    warn <- options(warn=-2)
    on.exit(options(warn)) # insure restore even in event of error
    nErrors <- 0
    N <- nrow(data)
    coefficients <- model$coeff
    coef.names <- names(coefficients)
    var.names <- model$var.names 
    ram <- model$ram 
    ram[coef.names, "start value"] <- coefficients 
    coefs <- matrix(numeric(0), R, length(coefficients))
    colnames(coefs) <- coef.names
    for (b in 1:R){
        for (try in 1:11){
            if (try > 10) stop("more than 10 consecutive convergence failures")
            res <- try(refit(), silent=TRUE)
            if (inherits(res, "try-error")) nErrors <- nErrors + 1
            else {
                coefs[b,] <- res
                break()
                }
            }
        }
    options(warn)
    if (nErrors > 0) warning("there were", nErrors, 
        "apparent convergence failures;\nthese are discarded from the",
        R, "bootstrap replications returned")
    res <- list(t0=coefficients, t=coefs, R=R, data=data, seed=seed,
        statistic=refit, sim="ordinary", stype="i", call=match.call(),
        strata=rep(1, N), weights=rep(1/N, N))
    class(res) <- c("bootsem", "boot")
    res
    }   
        
print.bootsem <- function(x, digits = getOption("digits"), ...){
    t <- x$t
    t0 <- x$t0
    result <- data.frame("Estimate"=t0, "Bias"=colMeans(t) - t0, 
        "Std.Error"=apply(t, 2, sd))
    rownames(result) <- names(t0)
    cat("Call: ")
    dput(x$call)
    cat("\n")
    print(result, digits=digits)
    invisible(x)
    }

summary.bootsem <- function(object,
    type=c("perc", "bca", "norm", "basic", "none"), level=0.95, ...){
    if ((!require("boot")) && (type != "none")) stop("boot package unavailable")
    type <- match.arg(type)
    t <- object$t
    t0 <- object$t0
    result <- data.frame("Estimate"=t0, "Bias"=colMeans(t) - t0, 
        "Std.Error"=apply(t, 2, sd))
    if (type != "none"){
        p <- length(t0)
        lower <- upper <- rep(0, p)
        low <- if (type == "norm") 2 else 4
        up  <- if (type == "norm") 3 else 5 
        for (i in 1:p){
            ci <- as.vector(boot.ci(object, type=type, index=i, 
                conf=level)[[type]])
            lower[i] <- ci[low]
            upper[i] <- ci[up]
            }
        result$Lower <- lower
        result$Upper <- upper
        }    
    rownames(result) <- names(t0)
    result <- list(table=result, call=object$call, level=level, type=type)
    class(result) <- "summary.bootsem"
    result
    }

print.summary.bootsem <- function(x, digits = getOption("digits"), ...){
    cat("Call: ")
    dput(x$call)
    cat("\n")
    if (x$type != "none") {
        cat(paste("Lower and upper limits are for the", 100*x$level, 
            "percent", x$type, "confidence interval\n\n"))
        }
    print(x$table, digits=digits)
    invisible(return(x))
    }
