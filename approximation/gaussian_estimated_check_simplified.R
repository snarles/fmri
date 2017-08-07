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

p <- dim(Omega)[1]

nms_s <- zeros(mc.reps, K)
y_s <- zeros(mc.reps, p)

for (i in 1:mc.reps) {
  mus <- randn(K, p)
  mu_hs <- mus + mvrnorm(K, mu = rep(0, p), Sigma = Xi)
  ys <- mus + mvrnorm(K, mu = rep(0, p), Sigma = Omega)
  y <- ys[1, ]
  delts <- t(t(mu_hs) - y)
  scores <- rowSums(delts^2)
  # record results
  nms_s[i, ] <- scores
  y_s[i, ] <- y
}

####
##  CHECK MEAN AND VARIANCE FOR SCORES Z_*, Z_1..., Z_{K-1}
####

aa <- -2 * p
bb <- 2 * TR2(Omega + Xi)
cc <- 2 * TR2(Omega)
dd <- 2 * TR2(2*eye(p) + Omega +  Xi)
ee <- 2 * TR2(eye(p) + Omega)
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

####
##  COMPARE ACTUAL MISCLASSIFICATION PROBS TO GAUSSIAN MISC PROBS
####

pmins <- do.call(pmin, data.frame(nms_s))
(mc_emp <- mean(nms_s[, 1] > pmins))
(mc_the <- 1 - pnormal_abcde(aa, bb, cc, dd, ee, K))
