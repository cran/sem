# last modified 2011-087-07 by J. Fox

read.moments <- function(...){
	.Deprecated("readMoments", package="sem")
	readMoments(...)
}

readMoments <- function(file="", diag=TRUE, 
        names=as.character(paste("X", 1:n, sep=""))){
    elements <- scan(file=file)
    m <- length(elements)
    d <- if (diag) 1 else -1
    n <- floor((sqrt(1 + 8*m) - d)/2)
    if (m != n*(n + d)/2) 
        stop("wrong number of elements (cannot make square matrix)")
    if (length(names) != n) stop("wrong number of variable names")
    X <- diag(n)
    X[upper.tri(X, diag=diag)] <- elements
    rownames(X) <- colnames(X) <- names
    t(X)
    }
    
