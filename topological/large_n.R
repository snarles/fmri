source("topological/large_n_source.R")
library(ade4)

####
## NUMERICAL CONVERGENCE TEST
####

n <- 1000; p <- 30; q <- 20
X <- randn(n, p)
G0 <- randn(p, q)
A0 <- randn(q, q); B0 <- randn(q, q)
sigma <- 1
truth <- list(G=G0, A=A0, B=B0)

# MA0 <- G0 %*% A0 %*% solve(t(A0) %*% A0, t(G0 %*% A0))
# MB0 <- G0 %*% B0 %*% solve(t(B0) %*% B0, t(G0 %*% B0))
# f2(MA0, MB0)
# f2(MA0, G0 %*% t(G0))

Y <- (X %*% G0  + sigma * randn(n, q)) %*% A0
W <- (X %*% G0  + sigma * randn(n, q)) %*% B0
sol <- jmle(X, Y, W, 100)
sol0 <- jmle(X, Y, W, 100, init=truth)
sep <- sep_mle_of(X, Y, W)

#f2(G0 %*% t(G0), sol0$G %*% t(sol0$G))
f2(G0 %*% t(G0), sol$G %*% t(sol$G))
#f2(G0 %*% t(G0))
f2(G0 %*% A0, sol$G %*% sol$A)
#f2(G0 %*% A0)
f2(G0 %*% B0, sol$G %*% sol$B)
#f2(G0 %*% B0)

f2(A0)
f2(sol$A)

truth %$% objf0(G, A, B)
(jof <- sol %$% objf0(G, A, B))
sol0 %$% objf0(G, A, B)

(lr <- 1/2*(jof - sep))

sol %$% log(det(t(A) %*% A))
truth %$% log(det(t(A) %*% A))

####
##  FUNCTIONS FOR COMPARISON TO MANTEL
####


simulate_mantel_jmle <- function(X, G0a, G0b = G0a) {
  n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(G0a)[2]
  A0 <- randn(q, q); B0 <- randn(q, q)
  sigma <- 1
  Y <- (X %*% G0a  + sigma * randn(n, q)) %*% A0
  W <- (X %*% G0b  + sigma * randn(n, q)) %*% B0
  dY <- dist(Y)
  dW <- dist(W)
  res <- mantel.randtest(dY, dW, 1e4)
  mantel_pv <- res$pvalue
  sol <- jmle(X, Y, W, 100)
  (jof <- sol %$% objf0(G, A, B))
  sep <- sep_mle_of(X, Y, W)
  (lr <- 1/2*(jof - sep))
  c(mantel_pv, lr)
}


####
##  Comparison with varying correlations
####

library(parallel)

## Fix X and G0a
n <- 1000; p <- 30; q <- 20
X <- randn(n, p)
G0a <- randn(p, q)
GE <- randn(p, q)
mc.its <- 7
mcc <- 7

## Independent case
G0b <- GE
res_ind <- do.call(rbind, mclapply(1:mc.its, 
              function(i) simulate_mantel_jmle(X, G0a, G0b), mc.cores = mcc))

## Identical case
G0b <- G0a
res_null <- do.call(rbind, mclapply(1:mc.its, 
              function(i) simulate_mantel_jmle(X, G0a, G0b), mc.cores = mcc))

## Correlated case
G0b <- sqrt(0.5) * G0a + sqrt(0.5) * GE
res_corr <- do.call(rbind, mclapply(1:mc.its, 
              function(i) simulate_mantel_jmle(X, G0a, G0b), mc.cores = mcc))

# MANTEL
boxplot.matrix(log(cbind(res_ind[, 1], res_corr[, 1], res_null[, 1])))

# JMLE
boxplot.matrix(log(cbind(res_ind[, 2], res_corr[, 2], res_null[, 2])))


