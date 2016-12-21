####
##  Precision paper
##  identification loss with low-dim, high-dim, and mixture
####

library(pracma)
library(lineId)

## low-dim

ntr <- 20
nte <- 100
p <- 3
q <- 5
sigma <- 0.5
B <- randn(p, q)

Xtr <- randn(ntr, p)
Ytr <- Xtr %*% B + sigma * randn(ntr, q)
Xte <- randn(nte, p)
Yte <- Xte %*% B + sigma * randn(nte, q)

Bhat <- solve(t(Xtr) %*% Xtr, t(Xtr) %*% Ytr)
Yhat <- Xte %*% Bhat
dmat <- pdist2(Yhat, Yte)
acs <- 1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
plot(acs)
lineId::piK()

