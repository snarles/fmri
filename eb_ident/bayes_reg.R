####
## Multivariate Bayesian Regression
####

library(magrittr)
library(pracma, warn.conflicts = FALSE)
library(MASS)

## Utility functions

f2 <- function(x) sum(x^2)

## Functions for Kronecker products

tkron_d_kron <- function(A, B, d) {
  a1 <- dim(A)[1]; a2 <- dim(A)[2]; b1 <- dim(B)[1]; b2 <- dim(B)[2]
  dmat <- matrix(d, b1, a1)
  # columns of dmat are diag(D1), diag(D2), ...
  C <- zeros(a2 * b2)
  if (a1 < b1) {
    bdbmat <- zeros(b2^2, a1)
    for (i in 1:a1) bdbmat[, i] <- as.numeric(t(B) %*% (dmat[, i] * B))
    for (i in 1:a2) {
      for (j in i:a2) {
        Cij <- matrix(bdbmat %*% (A[, i] * A[, j]), b2, b2)
        C[(i-1) * b2 + (1:b2), (j-1) * b2 + (1:b2)] <- Cij
        C[(j-1) * b2 + (1:b2), (i-1) * b2 + (1:b2)] <- t(Cij)
      }
    }
  } else {
    for (i in 1:a2) {
      for (j in i:a2) {
        dtemp <- as.numeric(dmat %*% (A[, i] * A[, j]))
        Cij <- t(B) %*% (dtemp * B)
        C[(i-1) * b2 + (1:b2), (j-1) * b2 + (1:b2)] <- Cij
        C[(j-1) * b2 + (1:b2), (i-1) * b2 + (1:b2)] <- t(Cij)
      }
    }  
  }
  C
}

kron_v <- function(A, B, cc) {
  a1 <- dim(A)[1]; a2 <- dim(A)[2]; b1 <- dim(B)[1]; b2 <- dim(B)[2]
  as.numeric(B %*% matrix(cc, b2, a2) %*% t(A))
}

## Posterior mean and covariance given a prior on Bvec
##  of the form Sigma_b %x% eye(pX)
##  and noise prior of form Sigma_e %x% Sigma_t

post_moments <- function(X, Y, Sigma_e, Sigma_b, Sigma_t = eye(dim(X)[1]), 
                         computeCov = TRUE, naive = FALSE, matrix = TRUE, ...) {
  n <- dim(X)[1]; pX <- dim(X)[2]; pY <- dim(Y)[2]
  Omega_e <- solve(Sigma_e)
  Omega_t <- solve(Sigma_t)
  Omega_b <- diag(1/diag(Sigma_b))
  yVec <- as.numeric(Y)
  xtx <- t(X) %*% Omega_t %*% X
  if (naive) {
    ans_mu <- solve(Omega_e %x% xtx + Omega_b %x% eye(pX),
                     t(eye(pY) %x% X) %*% (Omega_e %x% Omega_t) %*% yVec)
    if (matrix) ans_mu <- matrix(ans_mu, pX, pY)
    if (computeCov) {
      ans_cov <- solve(Omega_e %x% xtx + Omega_b %x% eye(pX))
      return(list(Mu = ans_mu, Cov = ans_cov))
    }
    return(ans_mu)
  }
  homega_b <- diag(1/sqrt(diag(Sigma_b)))
  hsigma_b <- diag(sqrt(diag(Sigma_b)))
  res <- eigen(hsigma_b %*% Omega_e %*% hsigma_b)
  D_e <- diag(res$values)
  V_e <- homega_b %*% res$vectors
  iV_e <- solve(V_e)
  resX <- eigen(xtx)
  V_x <- resX$vectors
  D_x <- diag(resX$values)
  d <- 1/(diag(D_e) %x% diag(D_x) + 1)
  temp <- kron_v(iV_e %*% Omega_e, t(V_x) %*% t(X) %*% Omega_t, yVec)
  ans_mu <- kron_v(t(iV_e), V_x, d * temp)
  if (matrix) ans_mu <- matrix(ans_mu, pX, pY)
  if (computeCov) {
    ans_cov <- tkron_d_kron(iV_e, t(V_x), d)
    return(list(Mu = ans_mu, Cov = ans_cov))
  }
  ans_mu
}


## Posterior predictive means
## Returns a list with entry for every row of Xte, which is a list with mu, cov

post_predictive <- function(X, Y, X_te, Sigma_e, Sigma_b, Sigma_t = eye(dim(X)[1]), 
                         naive = FALSE, ...) {
  n <- dim(X)[1]; pX <- dim(X)[2]; pY <- dim(Y)[2]
  L <- dim(X_te)[1]; ans <- as.list(numeric(L))
  Omega_e <- solve(Sigma_e)
  Omega_t <- solve(Sigma_t)
  Omega_b <- diag(1/diag(Sigma_b))
  yVec <- as.numeric(Y)
  xtx <- t(X) %*% Omega_t %*% X
  
  if (naive) {
    res <- post_moments(X, Y, Sigma_e, Sigma_b, Sigma_t, naive = TRUE, matrix = TRUE)
    B <- res$Mu
    SigmaB <- res$Cov
    for (i in 1:L) {
      x_star <- X_te[i, ]
      Ix <- eye(pY) %x% t(x_star)
      Mu <- as.numeric(t(x_star) %*% B)
      Cov <- Ix %*% SigmaB %*% t(Ix) + Sigma_e
      ans[[i]] <- list(Mu = Mu, Cov = Cov)
    }
  } else {
    homega_b <- diag(1/sqrt(diag(Sigma_b)))
    hsigma_b <- diag(sqrt(diag(Sigma_b)))
    res <- eigen(hsigma_b %*% Omega_e %*% hsigma_b)
    D_e <- diag(res$values)
    V_e <- homega_b %*% res$vectors
    iV_e <- solve(V_e)
    resX <- eigen(xtx)
    V_x <- resX$vectors
    D_x <- diag(resX$values)
    d <- 1/(diag(D_e) %x% diag(D_x) + 1)
    temp <- d * kron_v(iV_e %*% Omega_e, t(V_x) %*% t(X) %*% Omega_t, yVec)
    for (i in 1:L) {
      x_star <- X_te[i, ]
      Mu <- kron_v(t(iV_e), t(x_star) %*% V_x, temp)
      Cov <- tkron_d_kron(iV_e, t(V_x) %*% x_star, d) + Sigma_e
      ans[[i]] <- list(Mu = Mu, Cov = Cov)
    } 
  }
  ans
}