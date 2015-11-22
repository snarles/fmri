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
mce <- function(Omega, Xi = Omega, 
                A = solve(Id(Omega) + Omega - solve(Id(Omega) + Xi)),
                B = solve(Id(Omega) + Xi),
                K, mc.reps = 100) {
  p <- dim(Omega)[1]
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
mce2 <- function(Omega, Xi = Omega, 
                A = solve(Id(Omega) + Omega - solve(Id(Omega) + Xi)),
                B = solve(Id(Omega) + Xi),
                K, mc.reps = 100) {
  p <- dim(Omega)[1]
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

p <- 10
Omega <- 0.5 * cov(randn(2*p, p))
Xi <- 0.05 * cov(randn(2*p, p))
#A <- solve(Omega) + 0 * cov(randn(p, p))
A <- eye(p)
OmegaH <- cov(mvrnorm(100, mu = rep(0, p), Sigma = Omega))
XiH <- cov(mvrnorm(100, mu = rep(0, p), Sigma = Xi))
A <- solve(Id(Omega) + OmegaH - solve(Id(Omega) + XiH))
B <- solve(Id(Omega) + XiH)
K <- 100
mc.reps <- 1e4
1/2/sqrt(mc.reps)
mce(Omega, Xi, A, B, K, mc.reps=mc.reps)
mce2(Omega, Xi, A, B, K, mc.reps=mc.reps)

c(mce(Omega, Xi=0 * Omega, K=K), mc(solve(Omega), K))

