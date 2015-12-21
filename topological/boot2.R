library(pracma)
library(MASS)
f2 <- function(x, y = 0) sum((x-y)^2)

regression_data_model_ <- function(A, B, SigmaX, SigmaY) {
  p <- dim(A)[1]; q <- dim(A)[2]
  sampler <- function(nX, nY) {
    noise_mult <- 1
    if (nX == 0 && nY == 0) {
      nX <- 1; nY <- 1
      noise_mult <- 0
    }
    rawXc <- mvrnorm(nX, rep(0, q), SigmaX)
    rawXr <- rawXc %*% t(A) + randn(nX, p)
    rawYc <- mvrnorm(nY, rep(0, q), SigmaY)
    rawYr <- rawYc %*% t(B) + randn(nY, p)    
    dat <- rbind(cbind(0, rawXc, rawXr), cbind(1, rawYc, rawYr))
    list(p = p, q  = q, dat = dat)
  }
  # force eval
  lala <- sampler(10, 10)
  sampler
}

sample_moments <- function(res) {
  p <- res$p; q <- res$q; dat <- res$dat
  rawX <- dat[dat[,1] == 0, -1, drop = FALSE]
  rawY <- dat[dat[,1] == 1, -1, drop = FALSE]
  Xc <- rawX[, 1:q]
  Xr <- rawX[, -(1:q)]
  Yc <- rawY[, 1:q]
  Yr <- rawY[, -(1:q)]
  Ahat <- t(solve(t(Xc) %*% Xc, t(Xc) %*% Xr))
  Bhat <- t(solve(t(Yc) %*% Yc, t(Yc) %*% Yr))
  list(Ahat = Ahat, Bhat = Bhat)
}

p <- 5
q <- 2
SigmaX <- 3 * cov(randn(2*q, q)); SigmaY <- 3 * cov(randn(2 * q, q))
A_0 <- randn(p, q); B_0 <- svd(randn(p, p))$u %*% A_0
h0_small <- regression_data_model_(A_0, B_0, SigmaX, SigmaY)
dat <- h0_small(200, 200)
mus <- sample_moments(dat)
c(f2(mus$Ahat, A_0), f2(mus$Bhat, B_0))
