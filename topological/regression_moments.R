####
##  Moments of B' B
####

source('topological/rsa_boot_source.R')

library(magrittr)

gen_bhat <- function(n, B_true, SigmaX1, SigmaE1) {
  q <- dim(B_true)[1]; p <- dim(B_true)[2]
  X <- mvrnorm(n, rep(0, q), SigmaX1)
  Y <- X %*% B_true + mvrnorm(n, rep(0, p), SigmaE1)
  Bhat <- ginv(X) %*% Y
  Bhat  
}

####
##  Find the bias
####

p <- 5
q <- 2
B_true <- randn(q, p); A_0 <- t(B_true); B_0 <- svd(randn(p, p))$u %*% A_0
M_true <- B_true %*% t(B_true)
SigmaX1 <- 3 * cov(randn(2*q, q)); SigmaX2 <- 3 * cov(randn(2 * q, q))
SigmaE1 <- 2 * cov(randn(5 * p, p)); SigmaE2 <- cov(randn(5 * p, p))

n <- 200
mc.reps <- 10000
Bhat <- gen_bhat(n, B_true, SigmaX1, SigmaE1)
f2(Bhat, B_true)
f2(Bhat %*% t(Bhat), M_true)

Bhats <- lapply(1:mc.reps, function(i) gen_bhat(n, B_true, SigmaX1, SigmaE1))
BtB <- Bhats %>% lapply(function(v) v %*% t(v))
muM <- BtB %>% lapply(as.numeric) %>% do.call(what = rbind) %>% colMeans %>% matrix(q, q)
bias <- muM - M_true
the_bias <- solve(n * SigmaX1) * sum(diag(SigmaE1))
list(bias, the_bias)
