####
## factor model for generating probs
####

source("extrapolation/example2d.R")

library(pracma)

nX <- 200
nY <- 400
rk <- 3
xfax <- randn(nX, rk)
yfax <- randn(nY, rk)
ipmat <- xfax %*% t(yfax)
coefs <- c(-2, 1, -2, 3)
pmat <- exp(coefs[1] + coefs[2] * ipmat + coefs[3] * ipmat^2 + coefs[4] * ipmat)
for (i in 1:10) {
  rs <- rowSums(pmat)
  pmat <- pmat/rs
  pmat <- t(pmat)
}
avg_mc_acc_naive(pmat, 20)
avg_mc_acc_p(pmat, 20)
avg_mc_acc_naive(pmat, 40)
avg_mc_acc_p(pmat, 40)
rankconv <- (apply(pmat, 2, rank) - 0.5)/nrow(pmat)
plot(density(as.numeric(rankconv), weights = as.numeric(pmat)/sum(pmat)))
ps <- sample(as.numeric(rankconv), 10000, replace = TRUE, prob = as.numeric(pmat)/sum(pmat))
plot(sort(ps), type = "l")


