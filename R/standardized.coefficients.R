# last modified 25 July 2001 by J. Fox

standardized.coefficients <- function(object, digits=5){
    old.digits <- options(digits=digits)
    on.exit(options(old.digits))
    P <- object$P
    A <- object$A
    t <- object$t
    par <- object$coeff
    par.posn <- object$par.posn
    IAinv <- solve(diag(nrow(A)) - A)
    C <- IAinv %*% P %*% t(IAinv)
    ram <- object$ram
    par.names <- rep(' ', nrow(ram))
    for (i in 1:t) {
        which.par <- ram[,4] == i
        ram[which.par, 5] <- par[i]
        par.names[which.par] <- names(par)[i]
        }
    one.head <- ram[,1] == 1
    coeff <- ram[one.head, 5]
    coeff <- coeff *
        sqrt(diag(C[ram[one.head, 3], ram[one.head, 3]])/diag(C[ram[one.head, 2], ram[one.head, 2]]))
    var.names <- rownames(A)   
    par.code <- paste(var.names[ram[one.head,2]], c('<---', '<-->')[ram[one.head,1]],
                    var.names[ram[one.head,3]])
    coeff <- data.frame(par.names[one.head], coeff, par.code)
    names(coeff) <- c(" ", "Std. Estimate", " ")
    print(coeff, rowlab=rep(" ", nrow(coeff)))
    }

std.coef <- function (...){
    standardized.coefficients(...)
    }
