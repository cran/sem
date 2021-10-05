# The following functions were originally in the matrixcalc package
# written by Frederick Novomestky
# and distributed under the GPL license version 2 or higher

# last modified 2021-06-02

duplication.matrix <- function (n = 1) {
  return(D.matrix(n))
}

D.matrix <- function (n) {
  if (missing(n)) 
    stop("argument n is missing")
  if (!is.numeric(n)) 
    stop("argument n is not numeric")
  if (n != trunc(n)) 
    stop("argument n is not an integer")
  if (n < 2) 
    stop("argument n is less than 2")
  p <- n * (n + 1)/2
  nsq <- n * n
  Dt <- matrix(0, nrow = p, ncol = nsq)
  T <- T.matrices(n)
  u <- u.vectors(n)
  k <- u$k
  I <- u$I
  for (j in 1:n) {
    for (i in j:n) {
      Dt <- Dt + I[, k[i, j]] %*% t(vec(T[[i]][[j]]))
    }
  }
  return(t(Dt))
}

vech <- function (x) {
  if (!is.square.matrix(x)) 
    stop("argument x is not a square numeric matrix")
  return(t(t(x[!upper.tri(x)])))
}

is.square.matrix <- function (x) {
  if (!is.matrix(x)) 
    stop("argument x is not a matrix")
  return(nrow(x) == ncol(x))
}

T.matrices <- function (n) {
  if (missing(n)) 
    stop("argument n is missing")
  if (!is.numeric(n)) 
    stop("argument n is not numeric")
  if (n != trunc(n)) 
    stop("argument n is not an integer")
  if (n < 2) 
    stop("argument n is less than 2")
  E <- E.matrices(n)
  T <- list()
  for (i in 1:n) {
    T[[i]] <- list()
    for (j in 1:n) {
      if (i == j) {
        T[[i]][[j]] <- E[[i]][[j]]
      }
      else {
        T[[i]][[j]] <- E[[i]][[j]] + E[[j]][[i]]
      }
    }
  }
  return(T)
}

u.vectors <- function (n) {
  if (n != trunc(n)) 
    stop("argument n is not an integer")
  if (n < 2) 
    stop("argument n is less than 2")
  p <- n * (n + 1)/2
  I <- diag(rep(1, p))
  k <- matrix(0, nrow = n, ncol = n)
  for (j in 1:n) {
    for (i in j:n) {
      k[i, j] <- (j - 1) * n + i - 0.5 * j * (j - 1)
    }
  }
  return(list(k = k, I = I))
}

vec <- function (x) {
  if (!is.matrix(x)) {
    stop("argument x is not a matrix")
  }
  if (!is.numeric(x)) {
    stop("argument x is not a numeric matrix")
  }
  return(t(t(as.vector(x))))
}

E.matrices <- function (n) {
  if (missing(n)) 
    stop("argument n is missing")
  if (!is.numeric(n)) 
    stop("argument n is not numeric")
  if (n != trunc(n)) 
    stop("argument n is not an integer")
  if (n < 2) 
    stop("argument n is less than 2")
  I <- diag(rep(1, n))
  E <- list()
  for (i in 1:n) {
    E[[i]] <- list()
    for (j in 1:n) {
      E[[i]][[j]] <- I[i, ] %o% I[j, ]
    }
  }
  return(E)
}
