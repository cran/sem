# Two-Stage Least Squares
#   John Fox

# last modified 9 November 02 by J. Fox

tsls <- function(y, ...){
    UseMethod("tsls")
    }

tsls.default <- function (y, X, Z, names=NULL, ...) {
    n <- length(y)
    p <- ncol(X)
    invZtZ <- solve(crossprod(Z))
    XtZ <- crossprod(X, Z)
    V <- solve(XtZ %*% invZtZ %*% t(XtZ))
    b <- V %*% XtZ %*% invZtZ %*% crossprod(Z, y)
    residuals <- y - X %*% b
    s2 <- sum(residuals^2)/(n - p)
    V <- s2*V
    result<-list()
    result$n <- n
    result$p <- p
    b <- as.vector(b)
    names(b) <- names
    result$coefficients <- b
    rownames(V) <- colnames(V) <- names
    result$V <- V
    result$s <- sqrt(s2)
    result$residuals <- as.vector(residuals)
    result$response <- y
    result$model.matrix <- X
    result$instruments <- Z
    result
    }

    
tsls.formula <- function(formula, instruments, data, subset, 
    na.action, contrasts=NULL, ...){
    if (missing(na.action)) 
        na.action <- options()$na.action
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
        m$data <- as.data.frame(data)
    response.name <- deparse(formula[[2]])
    form <- as.formula(paste(response.name, '~', 
        deparse(formula[[3]]), '+', deparse(instruments[[2]])))
    m$formula <- form
    m$instruments <- m$formula <- m$contrasts <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, sys.frame(sys.parent()))
    na.act <- attr(mf, "na.action")
    Z <- model.matrix(instruments, data = mf, contrasts)
    response <- attr(attr(mf, "terms"), "response")
    y <- mf[,response]
    X <- model.matrix(formula, data=mf, contrasts)
    result <- tsls(y, X, Z, colnames(X))
    result$response.name <- response.name
    result$formula <- formula
    result$instruments <- instruments
    if (!is.null(na.act)) 
        result$na.action <- na.act
    class(result) <- "tsls"
    result
    }

print.tsls <- function(x, ...){
    cat("\nModel Formula: ")
    print(x$formula)
    cat("\nInstruments: ")
    print(x$instruments)
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\n")
    invisible(x)
    }
    
    
summary.tsls <- function(object, digits=4, ...){
    save.digits <- unlist(options(digits=digits))
    on.exit(options(digits=save.digits))
    cat("\n 2SLS Estimates\n")
    cat("\nModel Formula: ")
    print(object$formula)
    cat("\nInstruments: ")
    print(object$instruments)
    cat("\nResiduals:\n")
    print(summary(residuals(object)))
    cat("\n")
    df <- object$n - object$p
    std.errors <- sqrt(diag(object$V))
    b <- object$coefficients
    t <- b/std.errors
    p <- 2*(1 - pt(abs(t), df))
    table <- cbind(b, std.errors, t, p)
    rownames(table) <- names(b)
    colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)")
    print(table)
    cat(paste("\nResidual standard error:", round(object$s, digits),
        "on", df, "degrees of freedom\n\n"))
    }
    
residuals.tsls <- function(object, ...){
    res <- object$residuals
    if (is.null(object$na.action)) 
        res
    else naresid(object$na.action, res)
    }

coefficients.tsls <- function(object, ...){
    object$coefficients
    }
    
fitted.tsls <- function(object, ...){
    yhat <- as.vector(object$model.matrix %*% object$coefficients)
    if (is.null(object$na.action)) 
        yhat
    else napredict(object$na.action, yhat)
    }
    
anova.tsls <- function(object, model.2, s2, dfe, ...){
    if(class(model.2) != "tsls") stop('requires two models of class tsls')
    s2.1 <- object$s^2
    n.1 <- object$n 
    p.1 <- object$p
    dfe.1 <- n.1 - p.1
    s2.2 <- model.2$s^2
    n.2 <- model.2$n
    p.2 <- model.2$p
    dfe.2 <- n.2 - p.2
    SS.1 <- s2.1 * dfe.1
    SS.2 <- s2.2 * dfe.2
    SS <- abs(SS.1 - SS.2)
    Df <- abs(dfe.2 - dfe.1)
    if (missing(s2)){
        s2 <- if (dfe.1 > dfe.2) s2.1 else s2.2
        f <- (SS/Df) / s2
        RSS <- c(SS.1, SS.2)
        Res.Df <- c(dfe.1, dfe.2)
        SS <- c(NA, SS)
        P <- c(NA, 1 - pf(f, Df, min(dfe.1, dfe.2)))
        Df <- c(NA, Df)
        f <- c(NA, f)
        rows <- c("Model 1", "Model 2")
        }
    else{
        f <- (SS/Df) / s2
        RSS <- c(SS.1, SS.2, s2*dfe)
        Res.Df <- c(dfe.1, dfe.2, dfe)
        SS <- c(NA, SS, NA)
        P <- c(NA, 1 - pf(f, Df, min(dfe.1, dfe.2)), NA)
        Df <- c(NA, Df, NA)
        f <- c(NA, f, NA)
        rows <- c("Model 1", "Model 2", "Error")
        }
    table <- data.frame(Res.Df, RSS, Df, SS, f, P)
    head.1 <- paste("Model 1: ",format(object$formula), "  Instruments:", 
        format(object$instruments))
    head.2 <- paste("Model 2: ",format(model.2$formula), "  Instruments:", 
        format(model.2$instruments))
    names(table) <- c("Res.Df", "RSS", "Df", "Sum of Sq", "F", "Pr(>F)")
    row.names(table) <- rows
    structure(table, heading = c("Analysis of Variance", "", head.1, head.2, ""), 
        class = c("anova", "data.frame"))
    }
