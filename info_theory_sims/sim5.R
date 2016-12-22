####
##  Precision paper
##  identification loss with low-dim, high-dim, and mixture
####

library(pracma)
library(lineId)
source("info_theory_sims/fit_ident_curve.R")

## low-dim

ntr <- 20
nte <- 100
p <- 3
q <- 5
sigma <- 0.8
B <- randn(p, q)

Xtr <- randn(ntr, p)
Ytr <- Xtr %*% B + sigma * randn(ntr, q)
Xte <- randn(nte, p)
Yte <- Xte %*% B + sigma * randn(nte, q)

Bhat <- solve(t(Xtr) %*% Xtr, t(Xtr) %*% Ytr)
Yhat <- Xte %*% Bhat
dmat <- pdist2(Yhat, Yte)
acs <- 1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
acss <- sapply(1:30, function(i) {
  Xte <- randn(nte, p)
  Yte <- Xte %*% B + sigma * randn(nte, q)
  Yhat <- Xte %*% Bhat
  dmat <- pdist2(Yhat, Yte)
  1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
})
acs_mu <- rowMeans(acss)
plot(acs, type = "l")
k_sub <- 50
(I_implied <- fit_I_to_curve(acs[1:k_sub], wt_exp = 0))
(I_implied2 <- fit_I_to_curve(acs_mu[1:k_sub], wt_exp = 0, nits = 20, i_max = 4))

acs_hat <- acs_curve(I_implied, 1:nte)
acs_hat2 <- acs_curve(I_implied2, 1:nte)

plot(1:nte, acs, type = "l"); lines(acs_hat, col = "red"); abline(v = k_sub)
lines(acs_mu, col = "blue")
lines(acs_hat2, col = "green")



## high-dim

ntr <- 2000
nte <- 100
p <- 30
q <- 50
sigma <- 10
B <- randn(p, q)

Xtr <- randn(ntr, p)
Ytr <- Xtr %*% B + sigma * randn(ntr, q)
Xte <- randn(nte, p)
Yte <- Xte %*% B + sigma * randn(nte, q)

Bhat <- solve(t(Xtr) %*% Xtr, t(Xtr) %*% Ytr)
Yhat <- Xte %*% Bhat
dmat <- pdist2(Yhat, Yte)
acs <- 1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
acss <- sapply(1:30, function(i) {
  Xte <- randn(nte, p)
  Yte <- Xte %*% B + sigma * randn(nte, q)
  Yhat <- Xte %*% Bhat
  dmat <- pdist2(Yhat, Yte)
  1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
})
acs_mu <- rowMeans(acss)
plot(acs, type = "l")
k_sub <- 50
(I_implied <- fit_I_to_curve(acs[1:k_sub], wt_exp = 0))
(I_implied2 <- fit_I_to_curve(acs_mu[1:k_sub], wt_exp = 0))

acs_hat <- acs_curve(I_implied, 1:nte)
acs_hat2 <- acs_curve(I_implied2, 1:nte)

plot(1:nte, acs, type = "l"); lines(acs_hat, col = "red"); abline(v = k_sub)
lines(acs_mu, col = "blue")
lines(acs_hat2, col = "green")

## bi cluster

ntr <- 2000
nte <- 100
p <- 30
q <- 50
sigma <- 10
B <- randn(p, q)

Xtr <- randn(ntr, p)
Xtr[, 1] <- Xtr[, 1] + 15 * (rbinom(ntr, 1, .5) - .5)
Ytr <- Xtr %*% B + sigma * randn(ntr, q)
Xte <- randn(nte, p)
Xte[, 1] <- Xte[, 1] + 15 * (rbinom(nte, 1, .5) - .5)
Yte <- Xte %*% B + sigma * randn(nte, q)

Bhat <- solve(t(Xtr) %*% Xtr, t(Xtr) %*% Ytr)
Yhat <- Xte %*% Bhat
dmat <- pdist2(Yhat, Yte)
acs <- 1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
acss <- sapply(1:30, function(i) {
  Xte <- randn(nte, p)
  Yte <- Xte %*% B + sigma * randn(nte, q)
  Yhat <- Xte %*% Bhat
  dmat <- pdist2(Yhat, Yte)
  1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
})
acs_mu <- rowMeans(acss)
plot(acs, type = "l")
k_sub <- 50
(I_implied <- fit_I_to_curve(acs[1:k_sub], wt_exp = 0))
(I_implied2 <- fit_I_to_curve(acs_mu[1:k_sub], wt_exp = 0))

acs_hat <- acs_curve(I_implied, 1:nte)
acs_hat2 <- acs_curve(I_implied2, 1:nte)

plot(1:nte, acs, type = "l"); lines(acs_hat, col = "red"); abline(v = k_sub)
lines(acs_mu, col = "blue")
lines(acs_hat2, col = "green")




## predict 100th ac from 30, trying different wts
# k_sub <- 50
# 
# (I_unif <- fit_I_to_curve(acs[1:k_sub], wt_exp = 0))
# (ac_unif <- acs_curve(I_unif, 100))
# abs(ac_unif - acs[100])
# 
# (I_wt <- fit_I_to_curve(acs[1:k_sub], wt_exp = -1))
# (ac_wt <- acs_curve(I_wt, 100))
# abs(ac_wt - acs[100])

