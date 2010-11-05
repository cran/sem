# last modified 2011-11-04 by J. Fox


summary.objectiveML <- function(object, digits=5, conf.level=.90, robust=FALSE, analytic.se=TRUE, ...) {
	vcov <- vcov(object, robust=robust, analytic=analytic.se)
	if (any(is.na(vcov))) stop("coefficient covariances cannot be computed")
	norm.res <- normalizedResiduals(object)
	se <- sqrt(diag(vcov))
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
	chisqNull <- chisqNull(object)
	if(!robust) {
		chisq <- object$criterion * (N - (!object$raw))
	}
	else { 
		chisq <- object$adj.obj$chisq.scaled
		chisqNull <- chisqNull/object$adj.obj$c 
	}	
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
			if (is.na(res$objective) || res$objective < 0){
				max <- 0
				warning("cannot find upper bound of RMSEA")
				break
			}				
			if (sqrt(res$objective) < tail/100) break
			max <- max/2
		}
		lam.U <- if (max <= 1) NA else res$minimum
		max <- max(max, 1)
		while (max > 1){
			res <- optimize(function(lam) (1 - tail - pchisq(chisq, df, ncp=lam))^2, interval=c(0, max))
			if (sqrt(res$objective) < tail/100) break
			max <- max/2
			if (is.na(res$objective) || res$objective < 0){
				max <- 0
				warning("cannot find lower bound of RMSEA")
				break
			}				
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
		coeff <- data.frame(object$coeff, se, z, 2*pnorm(abs(z), lower.tail=FALSE), par.code)
		names(coeff) <- c("Estimate", "Std Error", "z value", "Pr(>|z|)", " ")
		row.names(coeff) <- names(object$coeff)
	}
	else coeff <- NULL
	AIC <- AIC(object)
	AICc <- AICc(object)
	BIC <- BIC(object)
	CAIC <- CAIC(object)
	SRMR <- sqrt(sum(standardizedResiduals(object)^2 * 
							upper.tri(diag(n), diag=TRUE))/(n*(n + 1)/2))
	ans <- list(chisq=chisq, df=df, chisqNull=chisqNull, dfNull=dfNull,
			GFI=GFI, AGFI=AGFI, RMSEA=RMSEA, NFI=NFI, NNFI=NNFI, CFI=CFI, BIC=BIC, SRMR=SRMR, 
			AIC=AIC, AICc=AICc, CAIC=CAIC, Rsq=Rsq(object),
			norm.res=norm.res, coeff=coeff, digits=digits, 
			iterations=object$iterations, aliased=object$aliased, raw=object$raw,
			robust=robust, robust.vcov=object$robust.vcov, adj.obj=object$adj.obj)
	class(ans) <- "summary.objectiveML"
	ans
}

print.summary.objectiveML <- function(x, ...){
	old.digits <- options(digits=x$digits)
	on.exit(options(old.digits))
	if (x$raw) cat("\nModel fit to raw moment matrix.\n")	
	if (x$robust && !is.null(x$robust.vcov)){
		cat("\n\nSatorra-Bentler Corrected Fit Statistics and Standard Errors:\n")
		cat("\n Adjusted Model Chisquare = ", x$adj.obj$chisq.scaled, "  Df = ", x$df, 
				"Pr(>Chisq) =", if (x$df > 0) pchisq(x$adj.obj$chisq.scaled, x$df, lower.tail=FALSE)
						else NA)		
		#use the scaled chisq for all other indices
		x$chisq<-x$adj.obj$chisq.scaled
		x$coeff <- x$coef
		x$coeff[,2] <- sqrt(diag(x$robust.vcov))
		x$coeff[,3] <- x$coeff[,1]/x$coeff[,2]
		x$coeff[,4] <- 2*pnorm(abs(x$coeff[,3]), lower.tail=FALSE)
		colnames(x$coeff)[2] <- "Corrected SE"
	}
	else{
		cat("\n Model Chisquare = ", x$chisq, "  Df = ", x$df, 
				"Pr(>Chisq) =", if (x$df > 0) pchisq(x$chisq, x$df, lower.tail=FALSE)
						else NA)
	}	
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
	if (!is.null(x$AIC)) cat("\n AIC = ", x$AIC)
	if (!is.null(x$AICc)) cat("\n AICc = ", x$AICc)
	if (!is.null(x$BIC)) cat("\n BIC = ", x$BIC)
	if (!is.null(x$CAIC)) cat("\n CAIC = ", x$CAIC, "\n")
	cat("\n Normalized Residuals\n")
	print(summary(as.vector(x$norm.res)))
	cat("\n R-square for Endogenous Variables\n")
	print(round(x$Rsq, 4))
	if (!is.null(x$coeff)){
		cat("\n Parameter Estimates\n")
		print(x$coeff, right=FALSE)
		if (!is.na(x$iterations)) cat("\n Iterations = ", x$iterations, "\n")
		if (!is.null(x$aliased)) cat("\n Aliased parameters:", x$aliased, "\n")
	}
	invisible(x)
}

summary.objectiveGLS <- function(object, ...){
	summary <- summary.objectiveML(object, ..., analytic.se=FALSE)
#	summary$chisqNull <- chisqNull(object) # object$chisqNull
	S <- object$S
	Sinv <- solve(S)
	C <- object$C
	SinvSmC <- Sinv %*% (S - C)
	SinvS <- Sinv %*% S
	n <- object$n
	summary$GFI <- 1 - sum(diag(SinvSmC %*% SinvSmC))/sum(diag(SinvS %*% SinvS))
	summary$AGFI <-  1 - (n*(n + 1)/(2*summary$df))*(1 - summary$GFI)
	summary
}

deviance.objectiveML <- function(object, ...) object$criterion * (object$N - (!object$raw))

df.residual.sem <- function(object, ...) {
	n.fix <- object$n.fix
	n <- object$n
	t <- object$t
	n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
}

Rsq <- function(model){
	A <- model$A
	P <- model$P
	IAinv <- solve(diag(nrow(A)) - A)
	C <- IAinv %*% P %*% t(IAinv)
	R2 <- 1 - diag(P)/diag(C)
	R2 <- R2[classifyVariables(model$semmod)$endogenous]
	R2
}
