# last modified 23 Dec 2001 by J. Fox

sem <- function(ram, ...){
    if (is.character(ram)) class(ram) <- 'mod'
    UseMethod('sem', ram)
    }

sem.mod <- function (ram, S, N, obs.variables=rownames(S), fixed.x=NULL, debug=F, ...){
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
    sem(ram=ram, S=S, N=N, param.names=pars, var.names=vars, fixed.x=fixed.x)
    }
     

sem.default <- function(ram, S, N, param.names=paste('Param', 1:t, sep=''), 
    var.names=paste('V', 1:m, sep=''), fixed.x=NULL, 
    analytic.gradient=T, heywood=F, warn=F, control=list()){
    is.triangular <- function(X) {
        is.matrix(X) && (nrow(X) == ncol(X)) && 
            (all(0 == X[upper.tri(X)])) || (all(0 == X[lower.tri(X)]))
        }    
    is.symmetric <- function(X) {
        is.matrix(X) && (nrow(X) == ncol(X)) && all(X == t(X))
        }
    if (is.triangular(S)) S <- S + t(S) - diag(diag(S))
    if (!is.symmetric(S)) stop('S must be a square triangular or symmetric matrix')
    con <- list(optim.control=list(),
        optim.method=if (heywood) "L-BFGS-B" else "BFGS", nlm.iterlim=100)
    if ((!is.matrix(ram)) | ncol(ram) != 5 | (!is.numeric(ram)))
        stop ('ram argument must be a 5-column numeric matrix')
    con[names(control)] <- control
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
    m <- max(ram[,2])
    t <- max(ram[,4])
    J <- matrix(0, n, m)
    J[cbind(1:n, observed)]<-1
    par.posn <- unlist(lapply(apply(outer(ram[,4], 1:t, '=='), 2, which), "[", 1))
    colnames(ram)<-c("heads", "to", "from", "parameter", "start value")
    rownames(ram)<-rep("",nrow(ram))
    rownames(ram)[par.posn]<-param.names
    fixed <- ram[,4] == 0
    sel.free <- ram[,4]
    sel.free[fixed] <- 1
    one.head <- ram[,1] == 1
    one.free <- which( (!fixed) & one.head )
    two.free <- which( (!fixed) & (!one.head) )
    arrows.1 <- ram[one.head, c(2,3)]
    arrows.2 <- ram[!one.head, c(2,3)]
    arrows.2t <- ram[!one.head, c(3,2)]
    arrows.1.free <- ram[one.free,c(2,3)]
    arrows.2.free <- ram[two.free,c(2,3)]
    sel.free.1 <- sel.free[one.free]
    sel.free.2 <- sel.free[two.free]
    start <- if (any(is.na(ram[,5][par.posn]))) startvalues(S, ram)
        else ram[,5][par.posn]
    bounds <- rep(-Inf, t)
    bounds[((ram[,1]==2) & (ram[,2]==ram[,3]))[par.posn]] <- if (heywood) .01 else -Inf
    objective.1 <- function(par){
        A <- P <- matrix(0, m, m)
        val <- ifelse (fixed, ram[,5], par[sel.free])
        A[arrows.1] <- val[one.head]
        P[arrows.2t] <- P[arrows.2] <- val[!one.head]
        I.Ainv <- solve(diag(m) - A)
        C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
        Cinv <- solve(C)
        F <- sum(diag(S %*% Cinv)) + log(det(C))
        F
        }
    objective.2 <- function(par){
        A <- P <- matrix(0, m, m)
        val <- ifelse (fixed, ram[,5], par[sel.free])
        A[arrows.1] <- val[one.head]
        P[arrows.2t] <- P[arrows.2] <- val[!one.head]
        I.Ainv <- solve(diag(m) - A)
        C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
        Cinv <- solve(C)
        F <- sum(diag(S %*% Cinv)) + log(det(C))
        grad.P <- t(I.Ainv) %*% t(J) %*% Cinv %*% (C - S) %*% Cinv %*% J %*% I.Ainv
        grad.A <- grad.P %*% P %*% t(I.Ainv)
        gradient <- rep(0, m)
        gradient[sel.free.1] <- grad.A[arrows.1.free]
        gradient[sel.free.2] <- grad.P[arrows.2.free]
        attributes(F) <- list(C=C, A=A, P=P, gradient=gradient)
        F
        }
    gradient <- function(par){
        A <- P <- matrix(0, m, m)
        val <- ifelse (fixed, ram[,5], par[sel.free])
        A[arrows.1] <- val[one.head]
        P[arrows.2t] <- P[arrows.2] <- val[!one.head]
        I.Ainv <- solve(diag(m) - A)
        C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)        
        Cinv <- solve(C)
        grad.P <- t(I.Ainv) %*% t(J) %*% Cinv %*% (C - S) %*% Cinv %*% J %*% I.Ainv
        grad.A <- grad.P %*% P %*% t(I.Ainv)
        gradient <- rep(0, m)
        gradient[sel.free.1] <- grad.A[arrows.1.free]
        gradient[sel.free.2] <- grad.P[arrows.2.free]
        gradient
        }
    if (!warn){
        save.warn <- options(warn=-1)
        on.exit(options(save.warn))
        }
    res <- if (con$optim.method == "L-BFGS-B") 
        optim(start, objective.1, method="L-BFGS-B",  
            gr=if (analytic.gradient) gradient,
            control=con$optim.control, lower=bounds)
        else optim(start, objective.1, method=con$optim.method,  
            gr=if (analytic.gradient) gradient,
            control=con$optim.control)
    convergence.1 <- res$convergence
    message.1 <- res$message
    coef.1 <- res$par
    if(res$convergence > 1) res$par <- start
    res <- nlm(if (analytic.gradient) objective.2 else objective.1, 
        res$par, hessian=T, iterlim=con$nlm.iterlim, check=F)
    convergence.2 <- res$code
    if (convergence.2 > 2) {
        coef.2 <- res$estimate
        res$par <- if (convergence.2 > 4) start else coef.2
        res <- nlm(if (!analytic.gradient) objective.2 else objective.1, 
            res$par, hessian=T, iterlim=con$nlm.iterlim, check=F)
        convergence.3 <- res$code
        }
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
    cov <- (2/(N - 1)) * solve(res$hessian)
    colnames(cov) <- rownames(cov) <- param.names
    result$cov <- cov
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
    result$convergence.1 <- convergence.1
    result$message.1 <- message.1
    result$coef.1 <- coef.1
    result$convergence.2 <- ultimate <- convergence.2
    if (exists("convergence.3", inherits=F)){
        result$coef.2 <- coef.2
        result$convergence.3 <- ultimate <- convergence.3
        }
    if (convergence.1 > 0) warning('initial optimization DID NOT converge')
    if (convergence.2 > 2) warning('second optimization DID NOT converge')
    if ((convergence.1 > 0) | (convergence.2 > 2)) 
        warning(paste('final optimization', if (ultimate > 2) 'DID NOT' else 'DID','converge'))
    class(result) <- "sem"
    result
    }
