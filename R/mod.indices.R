## last modified 3 Feb 2002 by J. Fox

mod.indices<-function(model, ...){
    UseMethod("mod.indices")
    }

mod.indices.sem <- function(model){        
    accumulate <- function(A, B, C, D, d) {
        res <- matrix(0, d^2, d^2)    
        for (ii in 1:(d^2)){
            for (jj in ii:(d^2)){
                g <- 1 + (ii - 1) %% d
                h <- 1 + (ii - 1) %/% d
                i <- 1 + (jj - 1) %% d
                j <- 1 + (jj - 1) %/% d
                res[ii, jj] <- res[jj, ii] <- A[g, i] * B[h, j] + C[g, j] * D[h, i]
                }
            }
        res
        }
    accumulate.asym <- function(A, B, C, D, d) {
        res <- matrix(0, d^2, d^2)    
        for (ii in 1:(d^2)){
            for (jj in 1:(d^2)){
                g <- 1 + (ii - 1) %% d
                h <- 1 + (ii - 1) %/% d
                i <- 1 + (jj - 1) %% d
                j <- 1 + (jj - 1) %/% d
                res[ii, jj] <- A[g, i] * B[h, j] + C[g, j] * D[h, i]
                }
            }
        res
        }        
    A <- model$A
    P <- model$P
    S <- model$S
    C <- model$C
    J <- model$J
    m <- model$m
    NM <- model$N - 1
    vars <- model$var.names    
    I.Ainv <- solve(diag(m) - A) 
    Cinv <- solve(C)    
    AA <- t(I.Ainv) %*% t(J)
    BB <- J %*% I.Ainv %*% P %*% t(I.Ainv)
    grad.P <- AA %*% Cinv %*% (C - S) %*% Cinv %*% t(AA)
    grad.A <- AA %*% Cinv %*% (C - S) %*% Cinv %*% BB 
    grad <- c(grad.A, grad.P)        
    dF.dBdB <- accumulate(AA %*% Cinv %*% t(AA), t(BB) %*% Cinv %*% BB,
                    AA %*% Cinv %*% BB, t(BB) %*% Cinv %*% t(AA), m)                
    dF.dPdP <- accumulate(AA %*% Cinv %*% t(AA), AA %*% Cinv %*% t(AA),
                    AA %*% Cinv %*% t(AA), AA %*% Cinv %*% t(AA), m)                
    dF.dBdP <- accumulate.asym(AA %*% Cinv %*% t(AA), t(BB) %*% Cinv %*% t(AA),
                    AA %*% Cinv %*% t(AA), t(BB) %*% Cinv %*% t(AA), m)    
    ram <- model$ram   
    fixed <- ram[,4] == 0
    one.head <- ram[,1] == 1
    one.free <- which( (!fixed) & one.head )
    two.free <- which( (!fixed) & (!one.head) )
    arrows.1.free <- ram[one.free,c(2,3)]
    arrows.2.free <- ram[two.free,c(2,3)]
    posn.matrix <- matrix(1:(m^2), m, m)
    posn.free <- c(posn.matrix[arrows.1.free], 
                (m^2) + posn.matrix[arrows.2.free])                    
    Hessian <- rbind( cbind(dF.dBdB,    dF.dBdP),
                      cbind(t(dF.dBdP), dF.dPdP))
    hessian <- Hessian[posn.free, posn.free]
    E.inv <- (2/NM) * solve(hessian)                      
    par.change <- mod.indices <- rep(0, 2*(m^2))                
    for (i in 1:(2*(m^2))) {
        k <- Hessian[i, i]
        d <- Hessian[i, posn.free]
     #   mod.indices[i] <- NM * grad[i]^2 / (k - d %*% E.inv %*% d)
        par.change[i] <- -grad[i]/ (k - d %*% E.inv %*% d)
        mod.indices[i] <- -NM * grad[i] * par.change[i]
        }
    mod.A <- matrix(mod.indices[1:(m^2)], m, m)
    mod.P <- matrix(mod.indices[-(1:(m^2))], m, m)
    par.A <- matrix(par.change[1:(m^2)], m, m)
    par.P <- matrix(par.change[-(1:(m^2))], m, m)
    mod.A[arrows.1.free] <- NA
    par.A[arrows.1.free] <- NA
    mod.P[arrows.2.free] <- NA
    par.P[arrows.2.free] <- NA
    rownames(mod.A) <- colnames(mod.A) <- vars
    rownames(mod.P) <- colnames(mod.P) <- vars
    rownames(par.A) <- colnames(par.A) <- vars
    rownames(par.P) <- colnames(par.P) <- vars
    result <- list(mod.A=mod.A, mod.P=mod.P, par.A=par.A, par.P=par.P)
    class(result) <- "sem.modind"
    result
    }

summary.sem.modind <- function(object, round=2, 
    print.matrices=c('both', 'par.change', 'mod.indices'), ...) {
    print.matrices <- match.arg(print.matrices)
    if (print.matrices != "mod.indices"){ 
        cat("\n Parameter change: A matrix\n")
        print(object$par.A)
        }
    if (print.matrices != "par.change"){ 
        cat("\n Modification indices: A matrix\n")
        print(round(object$mod.A, round))
        }
    if (print.matrices != "mod.indices"){ 
        cat("\n Parameter change: P matrix\n")
        print(object$par.P)
        }
    if (print.matrices != "par.change"){ 
        cat("\n Modification indices: P matrix\n")
        print(round(object$mod.P, round))
        }
    invisible(NULL)
    }

print.sem.modind <- function(x, n.largest=5, ...){
    cat("\n", n.largest, "largest modification indices, A matrix:\n")
    mod.A <- as.vector(x$mod.A)
    names <- rownames(x$mod.A)
    names(mod.A) <- outer(names, names, paste, sep=":")
    print(rev(sort(mod.A))[1:n.largest])
    cat("\n ", n.largest, "largest modification indices, P matrix:\n")
    mod.P <- as.vector(x$mod.P)
    names(mod.P) <- outer(names, names, paste, sep=":")
    print(rev(sort(mod.P[lower.tri(x$mod.P, diag=T)]))[1:n.largest])
    }
