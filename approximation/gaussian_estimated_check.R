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

mc.reps <- 1e4

cc <- 10
p <- 10; r <- 100; K <- 3
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

nms_s <- zeros(mc.reps, K)
y_s <- zeros(mc.reps, p)
Bmu_star_s <- zeros(mc.reps, p)
Bmu2_s <- zeros(mc.reps, p)

Ah <- sqrtm(A)
for (i in 1:mc.reps) {
  mus <- randn(K, p)
  mu_hs <- mus + mvrnorm(K, mu = rep(0, p), Sigma = Xi)
  ys <- mus + mvrnorm(K, mu = rep(0, p), Sigma = Omega)
  y <- ys[1, ]
  Bmu_hs <- mu_hs %*% t(B)
  delts <- t(t(Bmu_hs) - y)
  dA <- delts %*% Ah
  scores <- rowSums((delts %*% Ah)^2)
  scores2 <- scores - f2(Ah %*% y)
  # record results
  nms_s[i, ] <- scores
  y_s[i, ] <- y
  Bmu_star_s[i, ] <- Bmu_hs[1, ]
  Bmu2_s[i, ] <- Bmu_hs[2, ]
}

m_emp <- colMeans(nms_s)
m_the <- rep(TR(A %*% (eye(p) + Omega + B %*% (eye(p) + Omega/r) %*% B)), K)
m_the[1] <- m_the[1] - 2 * TR(A %*% B)

v_emp <- cov(nms_s)
c(v_emp[1,1],
  2 * TR2(A %*% (eye(p) + Omega + B %*% (eye(p) + Omega/r) %*% B - 2 * B)))
c(v_emp[1,2],
  2 * TR2(A %*% (eye(p) + Omega - B)))
c(v_emp[2,2],
  2 * TR2(A %*% (eye(p) + Omega + B %*% (eye(p) + Omega/r) %*% B)))
c(v_emp[2,3],
  2 * TR2(A %*% (eye(p) + Omega)))
