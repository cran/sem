# Last modified: 14 April 2007 by J. Fox

anova.sem <- function(object, model.2, ...){
    dev.1 <- deviance(object)
    df.1 <- df.residual(object)
    dev.2 <- deviance(model.2)
    df.2 <- df.residual(model.2)
    df <- abs(df.1 - df.2)
    if (df == 0) stop("the models have the same Df")
    if (object$N != model.2$N)
        stop("the models are fit to different numbers of observations")
    if ((nrow(object$S) != nrow(model.2$S)) || any(object$S != model.2$S))
        stop("the models are fit to different moment matrices")
    chisq <- abs(dev.1 - dev.2)
    table <- data.frame(c(df.1, df.2), c(dev.1, dev.2), c(NA, df), c(NA, chisq),
        c(NA, pchisq(chisq, df, lower.tail=FALSE)))
    names(table) <- c("Model Df", "Model Chisq", "Df", "LR Chisq", "Pr(>Chisq)")
    rownames(table) <- c("Model 1", "Model 2")
    structure(table, heading = c("LR Test for Difference Between Models", ""),
        class = c("anova", "data.frame"))
    }