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

mc.reps <- 1e4

cc <- 1
p <- 20; r <- 10; K <- 3
Omega <- 0.5 * cov(randn(10*p, p)) +  p * eye(p)
Omega <- Omega * TR(solve(Omega))/cc
#Omega <- eye(p) * p/cc
c(TR(solve(Omega)), TR2(solve(Omega)))
#Xi <- Omega/r
Xi <-  0.5 * cov(randn(10*p, p)) +  p * eye(p)
OmegaH <- cov(mvrnorm(K * r, rep(0, p), Omega)) 
#OmegaH <- Omega
XiH <- OmegaH/r
f2(Omega, OmegaH)

#A <- solve(eye(p) + OmegaH - solve(eye(p) + XiH))
A <- eye(p)
#B <- solve(eye(p) + XiH)
B <- eye(p)

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
bb <- 2 * TR2(A %*% (eye(p) + Omega + B %*% (eye(p) + Xi) %*% B - 2 * B))
cc <- 2 * TR2(A %*% (eye(p) + Omega - B))
dd <- 2 * TR2(A %*% (eye(p) + Omega + B %*% (eye(p) + Xi) %*% B))
ee <- 2 * TR2(A %*% (eye(p) + Omega))
m <- -aa/sqrt(dd- ee)
v <- (bb + ee - 2 * cc)/(dd - ee)

dim(nms_s)
emp_m <- colMeans(nms_s)
c(emp_m[1] - emp_m[2], aa)

Sigma <- diag(rep(dd - ee, K)) + ee
Sigma[1, 1] <- bb
Sigma[1, -1] <- cc; Sigma[-1, 1] <- cc
list(cov(nms_s), Sigma)

## formulas seem to work

gm <- eigen(A %*% B)$values
max(gm)
ls <- eigen(OmegaH)$values
c(aa, -2 * sum(1/(ls^2/r + (1 + 1/r) * ls)))
E <- Omega - OmegaH
c(bb, 2 * TR2(eye(p) + A %*% (E + B %*% E %*% B/r)))
c(cc, 2 * TR2(eye(p) + A %*% E))
c(dd, 2 * TR2(eye(p) + A %*% (E + 2 * B + B%*%E%*%B/r) ))
c(ee, 2 * TR2(eye(p) + A %*% (E + B) ))
c((bb + ee - 2 * cc)-(dd - ee), -4 * TR2(A %*% B) - 8/r*TR(A%*%B%*%A%*%B%*%E%*%B))


if (f2(Omega, OmegaH) == 0) {
  f2(A %*% (eye(p) + Omega), eye(p) + A %*% B)
  c(cc, bb,2 * p)
  c(ee, 2 * sum(((1 + ls) * (1 + ls/r)/((1 + ls)*(1 + ls/r) - 1))^2))
  c(dd, 2 * sum(( ( (1 + ls)*(1 + ls/r)+ 1 )/( (1 + ls)*(1+ls/r)- 1 )  )^2))
  c(ee, 2 * TR2(eye(p) + A %*% B))
  c(dd, 2 * TR2(eye(p) + 2*A %*% B))
  c(dd, 2 * sum((1 + 2 * gm)^2))
  c(m, 2*sum(gm)/sqrt(sum(4 * gm + 6 * gm^2)))
  c(v, 2 * sum(gm^2 + 2 * gm)/sum(6 * gm^2 + 4 * gm))
  # as gm -> 0
  c(m, sqrt(sum(gm)))
  c(v, 1)
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
