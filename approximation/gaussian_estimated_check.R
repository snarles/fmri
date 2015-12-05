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
source("approximation/pnorm_qnorm.R")

TR <- function(a) sum(diag(a))
TR2 <- function(a) sum(diag(a %*% a))

mc.reps <- 1e3

cc <- 1
p <- 10; r <- 100; K <- 3
Omega <- 0.5 * cov(randn(10*p, p)) + 1 * eye(p)
Omega <- Omega * TR(solve(Omega))/cc
#Omega <- eye(p) * p/cc
c(TR(solve(Omega)), TR2(solve(Omega)))
Xi <- Omega/r
#OmegaH <- cov(mvrnorm(K * r, rep(0, p), Omega)) 
OmegaH <- Omega
XiH <- OmegaH/r
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
  # record results
  nms_s[i, ] <- scores
  y_s[i, ] <- y
  Bmu_star_s[i, ] <- Bmu_hs[1, ]
  Bmu2_s[i, ] <- Bmu_hs[2, ]
}

####
##  CHECK MEAN AND VARIANCE FOR SCORES Z_*, Z_1..., Z_{K-1}
####

aa <- -2 * TR(A %*% B)
bb <- 2 * TR2(A %*% (eye(p) + Omega + B %*% (eye(p) + Omega/r) %*% B - 2 * B))
cc <- 2 * TR2(A %*% (eye(p) + Omega - B))
dd <- 2 * TR2(A %*% (eye(p) + Omega + B %*% (eye(p) + Omega/r) %*% B))
ee <- 2 * TR2(A %*% (eye(p) + Omega))

ls <- eigen(OmegaH)$values
c(aa, -2 * sum(1/(ls^2/r + (1 + 1/r) * ls)))
if (f2(Omega, OmegaH) == 0) {
  f2(A %*% (eye(p) + Omega), eye(p) + A %*% B)
  c(cc, bb,2 * p)
  c(ee, 2 * sum(((1 + ls) * (1 + ls/r)/((1 + ls)*(1 + ls/r) - 1))^2))
  c(dd, 2 * sum(( ( (1 + ls)*(1 + ls/r)+ 1 )/( (1 + ls)*(1+ls/r)- 1 )  )^2))
  c(ee, 2 * TR2(eye(p) + A %*% B))
  c(dd, 2 * TR2(eye(p) + 2*A %*% B))  
}

m_emp <- colMeans(nms_s)
m_the <- rep(TR(A %*% (eye(p) + Omega + B %*% (eye(p) + Omega/r) %*% B)), K)
m_the[1] <- m_the[1] + aa

v_emp <- cov(nms_s)
c(v_emp[1,1], bb)
c(v_emp[1,2], cc)
c(v_emp[2,2], dd)
c(v_emp[2,3], ee)

####
##  COMPARE ACTUAL MISCLASSIFICATION PROBS TO GAUSSIAN MISC PROBS
####

pmins <- do.call(pmin, data.frame(nms_s))
(mc_emp <- mean(nms_s[, 1] > pmins))
(mc_the <- 1 - pnormal_abcde(aa, bb, cc, dd, ee, K))
(mc_gauss <- 1 - pnormal_1min(m_emp, v_emp))
