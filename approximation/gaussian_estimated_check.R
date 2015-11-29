####
##  TESTING ASYMPTOTIC NORMALITY
####

library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc
library(lineId) ## use devtools::install('lineId')
source("approximation/gncx.R")

TR <- function(a) sum(diag(a))
TR2 <- function(a) sum(diag(a %*% a))

cc <- 10
p <- 100; r <- 2; K <- 1000
Omega <- 0.5 * cov(randn(20*p, p)) + 1000 * eye(p)
Omega <- Omega * TR(solve(Omega))/cc
c(TR(solve(Omega)), TR2(solve(Omega)))
Xi <- Omega/r
OmegaH <- cov(mvrnorm(K * r, rep(0, p), Omega)) 
XiH <- Omega/r
f2(Omega, OmegaH)

A <- solve(eye(p) + OmegaH - solve(eye(p) + XiH))
B <- solve(eye(p) + XiH)

p <- dim(Omega)[1]
Ah <- sqrtm(A)
mus <- randn(K, p)
mu_hs <- mus + mvrnorm(K, mu = rep(0, p), Sigma = Xi)
ys <- mus + mvrnorm(K, mu = rep(0, p), Sigma = Omega)
y <- ys[1, ]
Bmu_hs <- mu_hs %*% t(B)
delts <- t(t(Bmu_hs - y))
scores <- rowSums((delts %*% Ah)^2)
nms <- rowSums((Bmu_hs %*% Ah)^2)
c(mean(nms), TR(A %*% B %*% (eye(p) + Xi) %*% t(B)))
c(var(nms), 2 * TR2(A %*% B %*% (eye(p) + Xi) %*% t(B)))


dim(mu_hs)

lbls <- knn(mu_hs %*% t(B) %*% Ah, ys %*% Ah, cl=1:K)
sum(lbls != 1:K)/K

