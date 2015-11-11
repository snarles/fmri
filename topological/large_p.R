source("topological/large_n_source.R")
library(ade4)


####
##  FUNCTIONS FOR COMPARISON TO MANTEL
####


simulate_mantel_jmle <- function(X, G0a, G0b = G0a,
                                 sigma=1, lambda=1, alpha=0.5,
                                 Bsim = 0) {
  n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(G0a)[2]
  A0 <- randn(q, q); temp <- randn(q, q)
  Btemp <- pinv(G0b) %*% G0a %*% A0
  B0 <- sqrt(Bsim) * Btemp + sqrt(1-Bsim) * temp
  Y <- (X %*% G0a  + sigma * randn(n, q)) %*% A0
  W <- (X %*% G0b  + sigma * randn(n, q)) %*% B0
  dY <- dist(Y)
  dW <- dist(W)
  res <- mantel.randtest(dY, dW, 1e3)
  mantel_pv <- res$pvalue
  est <- init_est_reg(X, Y, W, lambda, alpha)
  #plot(as.numeric(est$M_Y), as.numeric(est$M_W))
  diff <- f2(est$M_Y, est$M_W)
  diff2 <- cor(as.numeric(est$M_Y), as.numeric(est$M_W))
  c(mantel_pv, diff, diff2)
}


####
##  Comparison with varying correlations
####

library(parallel)

## Fix X and G0a
n <- 200; p <- 500; q <- 200
sigma <- 0.01
Bsim <- 0.99
X <- randn(n, p)
G0a <- randn(p, q)
GE <- randn(p, q)
G0d <- G0a; G0d[, 1:q] <- 0
G0e <- GE; G0d[, -(1:q)] <- 0
mc.its <- 21
mcc <- 7

lambda <- 0.1; alpha <- 0.5
## Disjoint case
res_dis <- do.call(rbind, mclapply(1:mc.its, 
              function(i) 
                simulate_mantel_jmle(X, G0d, G0e,
                                    sigma, lambda, alpha, Bsim), mc.cores = mcc))

## Independent case
G0b <- GE
res_ind <- do.call(rbind, mclapply(1:mc.its, 
              function(i)
                simulate_mantel_jmle(X, G0a, G0b,
                                    sigma, lambda, alpha, Bsim), mc.cores = mcc))

## Identical case
G0b <- G0a
res_null <- do.call(rbind, mclapply(1:mc.its, 
              function(i) 
                simulate_mantel_jmle(X, G0a, G0b,
                                     sigma, lambda, alpha, Bsim), mc.cores = mcc))

## Correlated case
G0b <- sqrt(0.5) * G0a + sqrt(0.5) * GE
res_corr <- do.call(rbind, mclapply(1:mc.its, 
              function(i) 
                simulate_mantel_jmle(X, G0a, G0b,
                                    sigma, lambda, alpha, Bsim), mc.cores = mcc))

par(bg = "white")
# MANTEL
boxplot.matrix(log(cbind(
  dis=res_dis[,1], ind=res_ind[, 1], corr=res_corr[, 1], null=res_null[, 1])),
  main = "Mantel")

# diff
boxplot.matrix(log(cbind(
  dis=res_dis[,2], ind=res_ind[, 2], corr=res_corr[, 2], null=res_null[, 2])),
  main = "f2")

## diff2
boxplot.matrix(cbind(
  dis=res_dis[,3], ind=res_ind[, 3], corr=res_corr[, 3], null=res_null[, 3]),
  main = "cor")

