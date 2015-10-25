## Reference: Multivariate Statistics, by Fujikoshi Ulyanov and Shimizu

library(pracma)
library(MASS)

TR <- function(A) sum(diag(A))

eWAW <- function(n, Sigma, A = eye(dim(Sigma)[1]), empirical = 0) {
  p <- dim(Sigma)[1]
  if (empirical == 0) {
    return(
      n * Sigma %*% t(A) %*% Sigma + 
        n * TR(Sigma %*% A) * Sigma + 
        n^2 * Sigma %*% A %*% Sigma
    )
  } else {
    res <- array(0, dim = c(p, p, empirical))
    xall <- mvrnorm(n=n * empirical, mu=rep(0, p), Sigma = Sigma)
    for (i in 1:empirical) {
      xx <- xall[((i-1) * n) + (1:n), ]
      WW <- t(xx) %*% xx
      res[, , i] <- WW %*% A %*% WW
    }
    res2 <- apply(res, c(1, 2), mean)
    return(res2)
  }  
}


eTrAW <- function(n, Sigma, A, empirical = 0) {
  p <- dim(Sigma)[1]
  if (empirical == 0) {
    return(
      n * TR(A %*% Sigma)
    )
  }
  res <- numeric(empirical)
  xall <- mvrnorm(n=n * empirical, mu=rep(0, p), Sigma = Sigma)
  for (i in 1:empirical) {
    xx <- xall[((i-1) * n) + (1:n), ]
    WW <- t(xx) %*% xx
    res[i] <- TR(A %*% WW)
  }
  mean(res)
}

eTrAW_TrBW <- function(n, Sigma, A, W, empirical = 0){
  p <- dim(Sigma)[1]
  if (empirical == 0) {
    return(
      n * TR(A %*% Sigma %*% B %*% Sigma) +
        n * TR(t(A) %*% Sigma %*% B %*% Sigma) +
        n^2 * TR(A %*% Sigma) * TR(B %*% Sigma)
    )
  }
  res <- numeric(empirical)
  xall <- mvrnorm(n=n * empirical, mu=rep(0, p), Sigma = Sigma)
  for (i in 1:empirical) {
    xx <- xall[((i-1) * n) + (1:n), ]
    WW <- t(xx) %*% xx
    res[i] <- TR(A %*% WW) * TR(B %*% WW)
  }
  mean(res)
}

eTrAWBW <- function(n, Sigma, A, W, empirical = 0){
  p <- dim(Sigma)[1]
  if (empirical == 0) {
    return(
      n * TR(A %*% Sigma %*% t(B) %*% Sigma) +
        n * TR(A %*% Sigma) * TR(B %*% Sigma) +
        n^2 * TR(A %*% Sigma %*% B %*% Sigma)
    )
  }
  res <- numeric(empirical)
  xall <- mvrnorm(n=n * empirical, mu=rep(0, p), Sigma = Sigma)
  for (i in 1:empirical) {
    xx <- xall[((i-1) * n) + (1:n), ]
    WW <- t(xx) %*% xx
    res[i] <- TR(A %*% WW %*% B %*% WW)
  }
  mean(res)
}


eWiAWi <- function(n, Sigma, A, empirical = 0) {
  p <- dim(Sigma)[1]
  if (empirical == 0) {
    c2 <- 1/((n - p) * (n - p - 1) * (n - p - 3))
    c1 <- (n - p - 2) * c2
    ans <- c1 * solve(Sigma, A) %*% solve(Sigma) + 
      c2 * (solve(Sigma, t(A)) %*% solve(Sigma) +
              sum(diag(solve(Sigma, t(A)))) * solve(Sigma))
    return(ans)
  } else {
    res <- array(0, dim = c(p, p, empirical))
    xall <- mvrnorm(n=n * empirical, mu=rep(0, p), Sigma = Sigma)
    for (i in 1:empirical) {
      xx <- xall[((i-1) * n) + (1:n), ]
      WW <- t(xx) %*% xx
      res[, , i] <- solve(WW, A) %*% solve(WW)
    }
    res2 <- apply(res, c(1, 2), mean)
    return(res2)
  }
}

eTrWi2 <- function(n, Sigma) {
  p <- dim(Sigma)[1]
  c2 <- 1/((n - p) * (n - p - 1) * (n - p - 3))
  c1 <- (n - p - 2) * c2
  lambdas <- eigen(Sigma)$values
  ans <- c1 * sum(1/lambdas^2) + 
    c2 * (sum(1/lambdas^2) +sum(1/lambdas)^2)  
  ans
}


# "TESTS"
# n <- 10
# Sigma <- cov(randn(10, 5))
# A <- randn(5); B <- randn(5)
# list(eWAW(n, Sigma, A), eWAW(n, Sigma, A, empirical = 1e5))
# c(eTrAW(n, Sigma, A), eTrAW(n, Sigma, A, 1e5))
# c(eTrAW_TrBW(n, Sigma, A, B), eTrAW_TrBW(n, Sigma, A, B, 1e5))
# c(eTrAWBW(n, Sigma, A, B), eTrAWBW(n, Sigma, A, B, 1e5))
# list(eWiAWi(n, Sigma, A), eWiAWi(n, Sigma, A, 1e5))
# c(TR(eWiAWi(n, Sigma, eye(5))), eTrWi2(n, Sigma))

