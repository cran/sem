# last modified 9 Dec 04 by J. Fox

sem <- function(ram, ...){
    if (is.character(ram)) class(ram) <- 'mod'
    UseMethod('sem', ram)
    }

sem.mod <- function (ram, S, N, obs.variables=rownames(S), fixed.x=NULL, debug=FALSE, ...){
    parse.path <- function(path) {
        path.1 <- gsub('-', '', gsub(' ','', path))
        direction <- if (regexpr('<>', path.1) > 0) 2 
            else if (regexpr('<', path.1) > 0) -1
            else if (regexpr('>', path.1) > 0) 1
            else stop(paste('ill-formed path:', path))
        path.1 <- strsplit(path.1, '[<>]')[[1]]
        list(first=path.1[1], second=path.1[length(path.1)], direction=direction)
        }
    if ((!is.matrix(ram)) | ncol(ram) != 3) stop ('ram argument must be a 3-column matrix')
    startvalues <- as.numeric(ram[,3])
    par.names <- ram[,2]
    n.paths <- length(par.names)
    heads <- from <- to <- rep(0, n.paths)
    for (p in 1:n.paths){
        path <- parse.path(ram[p,1])
        heads[p] <- abs(path$direction)
        to[p] <- path$second
        from[p] <- path$first
        if (path$direction == -1) {
            to[p] <- path$first
            from[p] <- path$second
            }
        }
    ram <- matrix(0, p, 5)
    all.vars <- unique(c(to, from))
    latent.vars <- setdiff(all.vars, obs.variables)
    vars <- c(obs.variables, latent.vars)
    pars <- na.omit(unique(par.names))
    ram[,1] <- heads
    ram[,2] <- apply(outer(vars, to, '=='), 2, which)
    ram[,3] <- apply(outer(vars, from, '=='), 2, which)   
    par.nos <- apply(outer(pars, par.names, '=='), 2, which)
    ram[,4] <- unlist(lapply(par.nos, function(x) if (length(x) == 0) 0 else x))
    ram[,5]<- startvalues
    colnames(ram) <- c('heads', 'to', 'from', 'parameter', 'start')
    if (!is.null(fixed.x)) fixed.x <- apply(outer(vars, fixed.x, '=='), 2, which)
    n <- length(obs.variables)
    m <- length(all.vars)
    t <- length(pars)
    if (debug) {
        cat('\n observed variables:\n') 
        print(paste(paste(1:n,':', sep=''), obs.variables, sep=''))
        cat('\n')
        if (m > n){ 
            cat('\n latent variables:\n')
            print(paste(paste((n+1):m,':', sep=''), latent.vars, sep=''))
            cat('\n')
            }
        cat('\n parameters:\n') 
        print(paste(paste(1:t,':', sep=''), pars, sep=''))
        cat('\n\n RAM:\n')
        print(ram)
        }
    sem(ram=ram, S=S, N=N, param.names=pars, var.names=vars, fixed.x=fixed.x,
        debug=debug, ...)
    }
     

sem.default <- function(ram, S, N, param.names=paste('Param', 1:t, sep=''), 
    var.names=paste('V', 1:m, sep=''), fixed.x=NULL, debug=FALSE,
    analytic.gradient=TRUE, warn=FALSE, maxiter=500, par.size=c('ones', 'startvalues'), 
    refit=TRUE, start.tol=1E-6, ...){
    ord <- function(x) 1 + apply(outer(unique(x), x, "<"), 2, sum)
    is.triangular <- function(X) {
        is.matrix(X) && (nrow(X) == ncol(X)) && 
            (all(0 == X[upper.tri(X)])) || (all(0 == X[lower.tri(X)]))
        }    
    is.symmetric <- function(X) {
        is.matrix(X) && (nrow(X) == ncol(X)) && all(X == t(X))
        }
    if (is.triangular(S)) S <- S + t(S) - diag(diag(S))
    if (!is.symmetric(S)) stop('S must be a square triangular or symmetric matrix')
    if ((!is.matrix(ram)) | ncol(ram) != 5 | (!is.numeric(ram)))
        stop ('ram argument must be a 5-column numeric matrix')
    par.size <- if (missing(par.size)) {
        range <- range(diag(S))
        if (range[2]/range[1] > 100) 'startvalues' else 'ones'
        }
        else match.arg(par.size)
    n <- nrow(S)
    observed <- 1:n
    n.fix <- length(fixed.x)
    if (!is.null(fixed.x)){
        for (i in 1:n.fix){
            for (j in 1:i){
                ram <- rbind(ram, c(2, fixed.x[i], fixed.x[j], 
                    0, S[fixed.x[i], fixed.x[j]]))
                }
            }
        }
    m <- max(ram[,c(2,3)])
    missing.variances <- setdiff(1:m, ram[,2][ram[,2] == ram[,3]])
    if (length(missing.variances) > 0) warning(paste("The following variables have no variance or error-variance parameter (double-headed arrow):\n",
        paste(var.names[missing.variances], collapse=", "), 
        "\nThe model is almost surely misspecified; check also for missing covariances.\n"))
    t <- max(ram[,4])
    df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
    if (df < 0) stop(paste("The model has negative degrees of freedom =", df))
    J <- matrix(0, n, m)
    correct <- matrix(2, m, m)
    diag(correct) <- 1
    J[cbind(1:n, observed)]<-1
    par.posn <-  sapply(1:t, function(i) which(ram[,4] == i)[1])
    colnames(ram)<-c("heads", "to", "from", "parameter", "start value")
    rownames(ram)<-rep("",nrow(ram))
    rownames(ram)[par.posn]<-param.names
    fixed <- ram[,4] == 0
    sel.free <- ram[,4]
    sel.free[fixed] <- 1
    one.head <- ram[,1] == 1
    one.free <- which( (!fixed) & one.head )
    two.free <- which( (!fixed) & (!one.head) )
    arrows.1 <- ram[one.head, c(2,3), drop=FALSE]
    arrows.2 <- ram[!one.head, c(2,3), drop=FALSE]
    arrows.2t <- ram[!one.head, c(3,2), drop=FALSE]
    arrows.1.free <- ram[one.free,c(2,3), drop=FALSE]
    arrows.2.free <- ram[two.free,c(2,3), drop=FALSE]
    sel.free.1 <- sel.free[one.free]
    sel.free.2 <- sel.free[two.free]
    unique.free.1 <- unique(sel.free.1)
    unique.free.2 <- unique(sel.free.2)
    start <- if (any(is.na(ram[,5][par.posn]))) startvalues(S, ram, debug=debug, tol=start.tol)
        else ram[,5][par.posn]
    objective.1 <- function(par){
        A <- P <- matrix(0, m, m)
        val <- ifelse (fixed, ram[,5], par[sel.free])
        A[arrows.1] <- val[one.head]
        P[arrows.2t] <- P[arrows.2] <- val[!one.head]
        I.Ainv <- solve(diag(m) - A)
        C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
        Cinv <- solve(C)
        f <- sum(diag(S %*% Cinv)) + log(det(C))
        f
        }
    objective.2 <- function(par){
        A <- P <- matrix(0, m, m)
        val <- ifelse (fixed, ram[,5], par[sel.free])
        A[arrows.1] <- val[one.head]
        P[arrows.2t] <- P[arrows.2] <- val[!one.head]
        I.Ainv <- solve(diag(m) - A)
        C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
        Cinv <- solve(C)
        f <- sum(diag(S %*% Cinv)) + log(det(C))
        grad.P <- correct * t(I.Ainv) %*% t(J) %*% Cinv %*% (C - S) %*% Cinv %*% J %*% I.Ainv
        grad.A <- grad.P %*% P %*% t(I.Ainv)        
        gradient <- rep(0, t)
        gradient[unique.free.1] <- tapply(grad.A[arrows.1.free],ram[ram[,1]==1 & ram[,4]!=0, 4], sum)
        gradient[unique.free.2] <- tapply(grad.P[arrows.2.free],ram[ram[,1]==2 & ram[,4]!=0, 4], sum)
        attributes(f) <- list(C=C, A=A, P=P, gradient=gradient)
        f
        }
    if (!warn){
        save.warn <- options(warn=-1)
        on.exit(options(save.warn))
        }
    typsize <- if (par.size == 'startvalues') abs(start) else rep(1,t)
    res <- nlm(if (analytic.gradient) objective.2 else objective.1, 
        start, hessian=TRUE, iterlim=maxiter, print.level=if(debug) 2 else 0,
        typsize=typsize, ...)
    convergence <- res$code
    if (!warn) options(save.warn)
    par <- res$estimate
    names(par) <- param.names
    result <- list()
    result$var.names <- var.names
    obj <- objective.2(par)
    ram[par.posn, 5] <- start
    par.code <- paste(var.names[ram[,2]], c('<---', '<-->')[ram[,1]],
    var.names[ram[,3]])
    result$ram <- ram
    result$coeff <- par
    result$criterion <-  c(obj) - n - log(det(S))
    rownames(S) <- colnames(S) <- var.names[observed]
    result$S <- S
    result$J <- J
    C <- attr(obj, "C")
    rownames(C) <- colnames(C) <- var.names[observed]
    result$C <- C
    A <- attr(obj, "A")
    rownames(A) <- colnames(A) <- var.names
    result$A <- A
    P <- attr(obj, "P")
    rownames(P) <- colnames(P) <- var.names
    result$P <- P
    result$n.fix <- n.fix
    result$n <- n
    result$N <- N
    result$m <- m
    result$t <- t
    result$par.posn <- par.posn
    result$convergence <- convergence
    result$iterations <- res$iterations
    if (convergence > 2) 
        warning(paste('Optimization may not have converged; nlm return code = ',
            res$code, '. Consult ?nlm.\n', sep=""))
    qr.hess <- try(qr(res$hessian), silent=TRUE)
    if (class(qr.hess) == "try-error")
        warning("Could not compute QR decomposition of Hessian.\nOptimization probably did not converge.\n")
    else if (qr.hess$rank < t){ 
        warning(' singular Hessian: model is probably underidentified.\n')
        cov <- matrix(NA, t, t)
        colnames(cov) <- rownames(cov) <- param.names
        result$cov <- cov
        which.aliased <- qr.hess$pivot[-(1:qr.hess$rank)]
        aliased <- param.names[which.aliased]
        position.aliased <- is.element(ram[,4], which.aliased)
        if (refit){
            warning(' refitting without aliased parameters.\n')
            ram.refit <- ram[!position.aliased,]
            par.refit <- ram.refit[,4]
            par.refit[par.refit != 0] <- ord(par.refit[par.refit != 0])
            ram.refit[,4] <- par.refit
            iterations <- result$iterations
            result <- Recall(ram.refit, S, N, 
                param.names=param.names[-which.aliased], var.names=var.names,
                debug=debug, analytic.gradient=analytic.gradient, 
                warn=warn, maxiter=maxiter, par.size=par.size, refit=TRUE, ...)
            result$iterations <- iterations + result$iterations
            aliased <- c(aliased, result$aliased)
            }
        result$aliased <- aliased
        }
    else {
        cov <- (2/(N - 1)) * solve(res$hessian)
        colnames(cov) <- rownames(cov) <- param.names
        result$cov <- cov
        if (any(diag(cov) < 0)) {
            result$aliased <- c(param.names[diag(cov) < 0], result$aliased)
            warning("Negative parameter variances.\nModel is probably underidentified.\n")
            }
        }
    class(result) <- "sem"
    result
    }
