# last modified 25 July 2001 by J. Fox

sem <- function(S, ram, N, param.names=paste('Param', 1:t, sep=''), 
    var.names=paste('V', 1:m, sep=''), observed=1:n, fixed.x=NULL, 
    heywood=T, control=list()){
    con <- list(optim.control=list(), nlm.iterlim=100)
    con[names(control)] <- control
    n.fix <- length(fixed.x)
    if (!is.null(fixed.x)){
        for (i in 1:n.fix){
            for (j in 1:i){
                ram <- rbind(ram, c(2, fixed.x[i], fixed.x[j], 
                    0, S[fixed.x[i], fixed.x[j]]))
                }
            }
        }
    n <- nrow(S)
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
    start <- if (any(is.na(ram[,5][par.posn]))) start.values(S, ram, observed)
        else ram[,5][par.posn]
    bounds <- rep(-Inf, t)
    bounds[((ram[,1]==2) & (ram[,2]==ram[,3]))[par.posn]] <- if (heywood) .01 else -Inf
    objective <- function(par){
        A <- P <- matrix(0, m, m)
        val <- ifelse (fixed, ram[,5], par[sel.free])
        A[ram[one.head, c(2,3)]] <- val[one.head]
        P[ram[!one.head, c(2,3)]] <- P[ram[!one.head, c(3,2)]] <- val[!one.head]
        I.Ainv <- solve(diag(m) - A)
        C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
        F <- sum(diag(S %*% solve(C))) + log(det(C))
        attributes(F) <- list(C=C, A=A, P=P)
        F
        }
    res <- optim(start, objective, method="L-BFGS-B", 
        control=control$optim.control, lower=bounds)
    res <- nlm(objective, res$par, hessian=T, iterlim=con$nlm.iterlim)
    par <- res$estimate
    names(par) <- param.names
    result <- list()
    obj <- objective(par)
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
    result$convergence <- res$code
    if (result$convergence > 2) warning('did not converge')
    class(result) <- "sem"
    result
    }
