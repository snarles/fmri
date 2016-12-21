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
plot(acs, type = "l")

## predict 100th ac from 30, trying different wts

k_sub <- 50

(I_unif <- fit_I_to_curve(acs[1:k_sub], wt_exp = 0))
(ac_unif <- acs_curve(I_unif, 100))
abs(ac_unif - acs[100])

(I_wt <- fit_I_to_curve(acs[1:k_sub], wt_exp = -1))
(ac_wt <- acs_curve(I_wt, 100))
abs(ac_wt - acs[100])

