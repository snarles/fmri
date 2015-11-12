library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc
library(lineId) ## use devtools::install('lineId')
source("approximation/gncx.R")


TR <- function(a) sum(diag(a))

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

## uses conditioning on Y
mc2 <- function(Sigma, K, mc.reps = 1000) {
  p <- dim(Sigma)[1]
  Omega <- solve(Sigma)
  Amat <- eye(p) - solve(eye(p) + Omega)
  Ha_Y <- lineId::sqrtm(eye(p) + Omega)
  Ha_mu <- lineId::sqrtm(Amat)
  nmS <- function(v) t(v) %*% Sigma %*% v
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y <- as.numeric(Ha_Y %*% rnorm(p))
                  y_mu <- as.numeric(Ha_mu %*% rnorm(p) + Amat %*% y)
                  y_mus <- randn(p, K-1) + y
                  nms <- apply(y_mus, 2, nmS)
                  #nms <- colSums(Ha_S %*% y_mus)
                  min(nms) < t(y_mu) %*% Sigma %*% y_mu
                })
  mean(mcs)  
}

## uses conditioning on Y, draws random
mc2a <- function(Sigma, K, mc.reps = 1000) {
  p <- dim(Sigma)[1]
  Omega <- solve(Sigma)
  Amat <- eye(p) - solve(eye(p) + Omega)
  Ha_Y <- lineId::sqrtm(eye(p) + Omega)
  Ha_mu <- lineId::sqrtm(Amat)
  nmS <- function(v) t(v) %*% Sigma %*% v
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y <- as.numeric(Ha_Y %*% rnorm(p))
                  y_mu <- as.numeric(Ha_mu %*% rnorm(p) + Amat %*% y)
                  nms <- rgchisq0(K-1, Sigma, y)
                  min(nms) < t(y_mu) %*% Sigma %*% y_mu
                })
  mean(mcs)  
}

## use markov approximation
mc3 <- function(Sigma, K, mc.reps = 1000) {
  p <- dim(Sigma)[1]
  Omega <- solve(Sigma)
  Amat <- eye(p) - solve(eye(p) + Omega)
  Ha_Y <- lineId::sqrtm(eye(p) + Omega)
  Ha_mu <- lineId::sqrtm(Amat)
  nmS <- function(v) t(v) %*% Sigma %*% v
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y <- as.numeric(Ha_Y %*% rnorm(p))
                  y_mu <- as.numeric(Ha_mu %*% rnorm(p) + Amat %*% y)
                  mnms <- qlmb_gchisq(log(-log(1-runif(1))/(K-1)), Sigma, y)
                  mnms < t(y_mu) %*% Sigma %*% y_mu
                })
  mean(mcs)  
}

####
##  TESTS
####

p <- 10
Sigma <- cov(randn(2*p, p))
K <- 20
mc.reps <- 1e5
mc(Sigma, K, mc.reps)
mc2(Sigma, K, mc.reps)
mc2a(Sigma, K, mc.reps)
mc3(Sigma, K, 1000)

p <- 10
Sigma <- 5 * rexp(1) * cov(randn(2*p, p))
K <- 1e4
mc.reps <- 1e2
mc2(Sigma, K, mc.reps)
mc2a(Sigma, K, mc.reps)
mc3(Sigma, K, mc.reps)

