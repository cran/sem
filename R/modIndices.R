# last modified 2011-10-28 by J. Fox

mod.indices <- function(...){
	.Deprecated("modIndices", package="sem")
	modIndices(...)
}

modIndices <- function(model, ...){
    UseMethod("modIndices")
    }

# incorporates corrections by Michael Culbertson:

modIndices.objectiveML <- function(model, ...){        
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
	t <- model$t
	NM <- model$N - (!model$raw)
	vars <- model$var.names    
	I.Ainv <- solve(diag(m) - A) 
	Cinv <- solve(C)    
	AA <- t(I.Ainv) %*% t(J)
	BB <- J %*% I.Ainv %*% P %*% t(I.Ainv)
	correct <- matrix(2, m, m)
	diag(correct) <- 1
	grad.P <- correct * AA %*% Cinv %*% (C - S) %*% Cinv %*% t(AA)
	grad.A <- correct * AA %*% Cinv %*% (C - S) %*% Cinv %*% BB 
	grad <- c(grad.A, grad.P) * NM
	dF.dBdB <- accumulate(AA %*% Cinv %*% t(AA), t(BB) %*% Cinv %*% BB,
			AA %*% Cinv %*% BB, t(BB) %*% Cinv %*% t(AA), m)                
	dF.dPdP <- accumulate(AA %*% Cinv %*% t(AA), AA %*% Cinv %*% t(AA),
			AA %*% Cinv %*% t(AA), AA %*% Cinv %*% t(AA), m)                
	dF.dBdP <- accumulate.asym(AA %*% Cinv %*% t(AA), t(BB) %*% Cinv %*% t(AA),
			AA %*% Cinv %*% t(AA), t(BB) %*% Cinv %*% t(AA), m)    
	correct.BB <- correct.PP <- correct.BP <- matrix(1, m^2, m^2)
	d0 <- as.vector(diag(m) ==0 )
	d1 <- as.vector(diag(m) == 1)
	correct.BB[d0, d0] <- 2
	correct.PP[d1, d1] <- 0.5
	correct.PP[d0, d0] <- 2
	correct.BP[d0, d0] <- 2
	Hessian <- NM*rbind(cbind(dF.dBdB * correct.BB,    dF.dBdP * correct.BP),
			cbind(t(dF.dBdP * correct.BP), dF.dPdP * correct.PP))
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
	hessian <- Hessian[posn.free, posn.free]
	par.no <- ram[ram[, 4] > 0, 4]
	pars <- ram[, 4][!fixed]
	Z <- outer(1:t, pars, function(x, y) as.numeric(x == y))
	hessian <- Z %*% hessian %*% t(Z)
	E.inv <- solve(hessian)                      
	par.change <- mod.indices <- rep(0, 2*(m^2))       
	n.posn.free.1 <- length(posn.free) + 1
	for (i in 1:(2*(m^2))) {
		if (i %in% posn.free || qr(as.matrix(Hessian[c(posn.free, i), c(posn.free, i)]))$rank < n.posn.free.1){  
			par.change[i] <- mod.indices[i] <- NA  
		} 
		else {
			k <- Hessian[i, i]
			d <- Hessian[i, posn.free]
			d <- sapply(1:t, function(i) sum(d[which(par.no==i)]))
			par.change[i] <- -grad[i]/ (k - d %*% E.inv %*% d)
			mod.indices[i] <- -0.5 * grad[i] * par.change[i]
		}
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
	class(result) <- "modIndices"
	result
}


summary.modIndices <- function(object, round=2, 
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

print.modIndices <- function(x, n.largest=5, ...){
    cat("\n", n.largest, "largest modification indices, A matrix:\n")
    mod.A <- as.vector(x$mod.A)
    names <- rownames(x$mod.A)
    names(mod.A) <- outer(names, names, paste, sep="<-")
    print(rev(sort(mod.A))[1:n.largest])
    cat("\n ", n.largest, "largest modification indices, P matrix:\n")
    mod.P <- as.vector(x$mod.P)
    names(mod.P) <- outer(names, names, paste, sep="<->")
    print(rev(sort(mod.P[lower.tri(x$mod.P, diag=TRUE)]))[1:n.largest])
    }
