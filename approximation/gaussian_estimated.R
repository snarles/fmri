####
##  Misclassification rate with estimated centroids and covariance
####

library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc
library(lineId) ## use devtools::install('lineId')
source("approximation/gncx.R")
source("approximation/dgncx.R")

TR <- function(a) sum(diag(a))
Id <- function(a) eye(dim(a)[1])

## most naive implementation of mc()
mc <- function(Sigma, K, mc.reps = 1000) {
  p <- dim(Sigma)[1]
  mcs <- sapply(1:mc.reps,
                function(i) {
                  mus <- mvrnorm(n = K, mu = rep(0, p), Sigma = Sigma)
                  ys <- mus + randn(K, p)
                  lbls <- knn(mus, ys, cl=1:K)
                  sum(lbls != 1:K)/K
                })
  mean(mcs)
}

## most naive implementation of mc()
mce <- function
(Sigma, Omega = Id(Sigma), Xi = Omega,
 SigmaH = Sigma, OmegaH = Omega, XiH = Xi,
 K, mc.reps = 100) {
  p <- dim(Sigma)[1]
  A <- solve(SigmaH + OmegaH - SigmaH %*% solve(SigmaH + XiH) %*% SigmaH)
  B <- SigmaH %*% solve(SigmaH + XiH)
  Ah <- sqrtm(A)
  mcs <- sapply(1:mc.reps,
                function(i) {
                  mus <- mvrnorm(K, mu = rep(0, p), Sigma = Sigma)
                  mu_hs <- mus + mvrnorm(K, mu = rep(0, p), Sigma = Xi)
                  ys <- mus + mvrnorm(K, mu = rep(0, p), Sigma = Omega)
                  lbls <- knn(mu_hs %*% t(B) %*% Ah, ys %*% Ah, cl=1:K)
                  sum(lbls != 1:K)/K
                })
  mean(mcs)
}

## Converts to Sigma = I
mce1 <- function
(Sigma, Omega = Id(Sigma), Xi = Omega,
 SigmaH = Sigma, OmegaH = Omega, XiH = Xi,
 K, mc.reps = 100) {
  p <- dim(Omega)[1]
  ## convert everything if needed
  if (f2(Sigma, eye(p)) > 1e-8) {
    si <- isqrtm(Sigma)
    Omega <- si %*% Omega %*% si
    Xi <- si %*% Xi %*% si
    SigmaH <- si %*% SigmaH %*% si
    XiH <- si %*% Xi %*% si    
  }
  ## algo with Sigma = I
  A <- solve(SigmaH + OmegaH - SigmaH %*% solve(SigmaH + XiH) %*% SigmaH)
  B <- SigmaH %*% solve(SigmaH + XiH)
  Ah <- sqrtm(A)
  mcs <- sapply(1:mc.reps,
                function(i) {
                  mus <- randn(K, p)
                  mu_hs <- mus + mvrnorm(K, mu = rep(0, p), Sigma = Xi)
                  ys <- mus + mvrnorm(K, mu = rep(0, p), Sigma = Omega)
                  lbls <- knn(mu_hs %*% t(B) %*% Ah, ys %*% Ah, cl=1:K)
                  sum(lbls != 1:K)/K
                })
  mean(mcs)
}

## uses conditioning on Y
mce2 <- function
(Sigma, Omega = Id(Sigma), Xi = Omega,
 SigmaH = Sigma, OmegaH = Omega, XiH = Xi,
 K, mc.reps = 100) {
  p <- dim(Omega)[1]
  ## convert everything if needed
  if (f2(Sigma, eye(p)) > 1e-8) {
    si <- isqrtm(Sigma)
    Omega <- si %*% Omega %*% si
    Xi <- si %*% Xi %*% si
    SigmaH <- si %*% SigmaH %*% si
    XiH <- si %*% Xi %*% si    
  }
  ## algo with Sigma = I
  A <- solve(SigmaH + OmegaH - SigmaH %*% solve(SigmaH + XiH) %*% SigmaH)
  B <- SigmaH %*% solve(SigmaH + XiH)
  Ah <- sqrtm(A)
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y <- mvrnorm(1, rep(0, p), Sigma = eye(p) + Omega)
                  muh_star <- solve(eye(p) + Omega, y) +
                    mvrnorm(1, rep(0, p), eye(p) + Xi - solve(eye(p) + Omega))
                  muhs <- mvrnorm(K-1, rep(0, p), eye(p) + Xi)
                  muhsB <- B %*% cbind(muh_star, t(muhs))
                  diffs <- colSums((Ah %*% (muhsB - y))^2)
                  order(diffs)[1]!=1
                })
  mean(mcs)
}

## uses conditioning on Y, draws random
mce2a <- function
(Sigma, Omega = Id(Sigma), Xi = Omega,
 SigmaH = Sigma, OmegaH = Omega, XiH = Xi,
 K, mc.reps = 100) {
  p <- dim(Omega)[1]
  ## convert everything if needed
  if (f2(Sigma, eye(p)) > 1e-8) {
    si <- isqrtm(Sigma)
    Omega <- si %*% Omega %*% si
    Xi <- si %*% Xi %*% si
    SigmaH <- si %*% SigmaH %*% si
    XiH <- si %*% Xi %*% si    
  }
  ## algo with Sigma = I
  A <- solve(SigmaH + OmegaH - SigmaH %*% solve(SigmaH + XiH) %*% SigmaH)
  B <- SigmaH %*% solve(SigmaH + XiH)
  Ah <- sqrtm(A)
  cov_star <- B %*% (eye(p) + Xi - solve(eye(p) + Omega)) %*% t(B)
  cov_i <- B %*% (eye(p) + Xi) %*% t(B)
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y <- mvrnorm(1, rep(0, p), Sigma = eye(p) + Omega)
                  mu_star <- ((eye(p) - B %*% solve(eye(p) + Omega)) %*% y)[,1]
                  nm_star <- rdgchisq0(1, cov_star, mu_star, A)
                  nms <- rdgchisq0(K-1, cov_i, y, A)
                  min(nms) < nm_star
                })
  mean(mcs)
}

####
##  Basic tests
####

p <- 10
Sigma0 <- cov(randn(5 *p, p))
Omega0 <- 0.5 * cov(randn(2*p, p))
Xi0 <- 0.05 * cov(randn(2*p, p))

Sigma <- eye(p); Omega <- Omega0; Xi <- 0 * Omega
SigmaH <- Sigma; OmegaH <- Omega; XiH <- Xi

## check equivalences between mc and mce
K <- 200
c(mc(solve(Omega0), K),
  mce(eye(p), Omega0, Xi=0 * Omega0, K=K), 
  mce(solve(Omega0), eye(p), Xi=0 * Omega0, K=K),
  mce1(eye(p), Omega0, Xi=0 * Omega0, K=K), 
  mce1(solve(Omega0), eye(p), Xi=0 * Omega0, K=K))

## check mce, mce1, mce2
nest <- 20
SigmaH <- cov(mvrnorm(nest, mu = rep(0, p), Sigma = Sigma))
OmegaH <- cov(mvrnorm(nest, mu = rep(0, p), Sigma = Omega))
XiH <- cov(mvrnorm(nest, mu = rep(0, p), Sigma = Xi))
K <- 100
mc.reps <- 1e4
1/2/sqrt(mc.reps)
mce(Sigma, Omega, Xi, SigmaH, OmegaH, XiH, K, mc.reps=mc.reps)
mce1(Sigma, Omega, Xi, SigmaH, OmegaH, XiH, K, mc.reps=mc.reps)
mce2(Sigma, Omega, Xi, SigmaH, OmegaH, XiH, K, mc.reps=mc.reps)
mce2a(Sigma, Omega, Xi, SigmaH, OmegaH, XiH, K, mc.reps=mc.reps)


