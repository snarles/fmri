####
## Bayesian binomial estimation
####


## frequentist unbiased
binmom <- function(succ, tot, k) {
  choose(succ, k)/choose(tot, k)
}

## Bayesian
bbinom <- function(succ, tot, k, a = 1, b = 1) {
  a2 = a + succ; b2 = b + tot - succ
  beta(a2 + k, b2)/beta(a2, b2)
}


## test it
# a <- 2
# b <- 3
# ps <- rbeta(10000, a, b)
# ps <- sort(ps)
# c(mean(ps), a/(a + b))
# 
# k <- 10
# Ys <- rbinom(length(ps), k, prob = ps)
# K <- 7
# ests1 <- sapply(Ys, function(v) binmom(v, k, K))
# plot(ps^K, ests1, pch = ".")
# psK <- ps^K
# lm(ests1 ~ psK)
# ests2 <- sapply(Ys, function(v) bbinom(v, k, K, a, b))
# plot(ps^K, ests2, pch = ".")
# sum((ests1 - psK)^2)
# sum((ests2 - psK)^2)
# 
# K <- 20
# psK <- ps^K
# ests2 <- sapply(Ys, function(v) bbinom(v, k, K, a, b))
# sum((ests2 - ps^K)^2)
# res <- ksmooth(psK, ests2, bandwidth = 0.1)
# plot(res, type = "l")
# points(ps^K, ests2, pch = ".")
# abline(0, 1)
