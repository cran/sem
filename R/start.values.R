# last modified 14 April 2001 by J. Fox

start.values <- function(S, ram, observed){
    n <- nrow(S)        
    m <- max(ram[,2])            
    t <- max(ram[,4])   
    s <- sqrt(diag(S))
    R <- S/outer(s,s)
    latent<-(1:m)[-observed]
    par.posn <- unlist(lapply(apply(outer(ram[,4], 1:t, '=='), 2, which), "[", 1))
    one.head <- ram[,1] == 1
    start <- (ram[,5])[par.posn]
    A.pat <-matrix(FALSE, m, m)
    A.pat[ram[one.head, c(2,3)]] <- TRUE
    P.pat <- C <- matrix(0, m, m)
    P.pat[ram[!one.head, c(2,3)]] <- P.pat[ram[!one.head, c(3,2)]] <- 1
    C[observed, observed] <- R
    for (l in latent) {
        indicators <- A.pat[observed, l]
        for (j in observed){
            C[j, l] <- C[l, j] <- if (!any(indicators)) runif(1, .3, .5)
                else {   
                        numerator <- sum(R[j, observed[indicators]])
                        denominator <- sqrt(sum(R[observed[indicators], observed[indicators]]))
                        numerator/denominator
                    }
            }
        }
    for (l in latent){
        for (k in latent){
            C[l, k] <- if (l==k) 1 else {
                                indicators.l <- A.pat[observed, l]
                                indicators.k <- A.pat[observed, k]
                                if ((!any(indicators.l)) | (!any(indicators.k))) runif(1, .3, .5) else {
                                    numerator <- sum(R[observed[indicators.l], observed[indicators.k]])
                                    denominator <- sqrt( sum(R[observed[indicators.l], observed[indicators.l]])
                                        * sum(R[observed[indicators.k], observed[indicators.k]]))
                                    numerator/denominator}
                                    }
            }
        }
    A <- matrix(0, m, m)
    for (j in 1:m){
        ind <- A.pat[j,]
        if (!any(ind)) next
        A[j, ind] <- inv(C[ind, ind]) %*% C[ind, j]
        }
    A[observed,] <- A[observed,]*matrix(s, n, m)
    A[,observed] <- A[,observed]*matrix(s, m, n, byrow=T)
    C[observed,] <- C[observed,]*matrix(s, n, m)
    C[,observed] <- C[,observed]*matrix(s, m, n, byrow=T)
    P <- (diag(m) - A) %*% C %*% t(diag(m) - A)
    P <- P.pat * P
    for (par in 1:t){
        if (!is.na(start[par])) next
        posn <- par.posn[par]
        if (ram[posn, 1] == 1) start[par] <- A[ram[posn, 2], ram[posn, 3]]
            else start[par] <- P[ram[posn, 2], ram[posn, 3]]
        }
    start
    }
