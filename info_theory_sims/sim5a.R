####
##  Precision paper
##  get a 'good example' for identification loss
####

library(pracma)
library(lineId)
source("info_theory_sims/fit_ident_curve.R")
#source("idloss/idLoss.R")
library(glmnet)
library(parallel)
mcc <- 3
## sparsity

set.seed(0)

ntr <- 1000
nte <- 1000
p <- 200
q <- 200
sigma <- 3
B <- randn(p, q) * matrix(rbinom(p * q, 1, 0.02), p, q)

Xtr <- randn(ntr, p)
Ytr <- Xtr %*% B + sigma * randn(ntr, q)
Xte <- randn(nte, p)
Yte <- Xte %*% B + sigma * randn(nte, q)

t1 <- proc.time()
enet_res <- mclapply(1:ncol(Ytr), function(i) {
  cv.glmnet(Xtr, Ytr[, i], alpha = 0)
}, mc.cores = mcc)
proc.time() - t1

Yhat <- sapply(enet_res, function(r) predict(r, Xte))

t2 <- proc.time()
dmat <- pdist2(Yhat, Yte)
proc.time()- t2

acs <- 1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
acs

acss <- mclapply(1:30, function(i) {
  Xte <- randn(nte, p)
  Yte <- Xte %*% B + sigma * randn(nte, q)
  #Yhat <- Xte %*% Bhat
  Yhat <- sapply(enet_res, function(r) predict(r, Xte))
  dmat <- pdist2(Yhat, Yte)
  1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
}, mc.cores = mcc)
acss <- do.call(cbind, acss)
acs_mu <- rowMeans(acss)
#plot(acs, type = "l")
k_sub <- 50
(I_implied <- fit_I_to_curve(acs[1:k_sub], wt_exp = 0))
(I_implied2 <- fit_I_to_curve(acs_mu[1:k_sub], wt_exp = 0, nits = 20, i_max = 10))

acs_hat <- acs_curve(I_implied, 1:nte)
acs_hat2 <- acs_curve(I_implied2, 1:nte)

plot(1:nte, acs, type = "l", ylim = c(0, 1),
     xlab = "k", ylab = "accuracy")
lines(acs_hat, col = "red"); abline(v = k_sub)
lines(acs_mu, col = "blue")
#lines(acs_hat2, col = "green")

m_err <- sqrt(mean(((acs_hat - acs)^2)[-(1:k_sub)]))
m_err
m_err2 <- sqrt(mean(((acs_hat - acs_mu)^2)[-(1:k_sub)]))
m_err2
#title(sub = paste("RMSE =", m_err2))
