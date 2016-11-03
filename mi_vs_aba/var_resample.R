###
##  Variance of average GA when resampling K choose k classes
###

library(pracma)
library(lineId)
bigK <- 1000
K <- 50
k <- 3

plikes_pop <- randn(bigK) + 2 * eye(bigK)

resample_accs <- function(plikes, k, mc.reps = 100) {
  K <- nrow(plikes)
  sapply(1:mc.reps, function(i) {
    ksamp <- sample(K, k)
    psub <- plikes[ksamp, ksamp]
    maxr <- apply(psub, 1, max)
    sum(diag(psub) == maxr)/k
  })
}

accscore <- function(pl) {
  sum(diag(pl) == apply(pl, 1, max))/nrow(pl)
}

cov_accs_overlap <- function(plikes, k, ol = 1, mc.reps = 100) {
  K <- nrow(plikes)
  res <- t(sapply(1:mc.reps, function(i) {
    ksamp0 <- sample(K, ol)
    ksamp1 <- sample(setdiff(1:K, ksamp0), k - ol)
    ksamp2 <- sample(setdiff(1:K, c(ksamp0, ksamp1)), k - ol)
    ksampA <- c(ksamp0, ksamp1)
    ksampB <- c(ksamp0, ksamp2)
    psubA <- accscore(plikes[ksampA, ksampA])
    psubB <- accscore(plikes[ksampB, ksampB])
    c(psubA, psubB)
  }))
  dim(res)
  mean(res)
  cov(res)[1, 2]
}


## actual variance
mc.reps <- 1e4
aba_hats <- sapply(1:mc.reps, function(i) {
  Ksamp <- sample(bigK, K)
  plikes <- plikes_pop[Ksamp, Ksamp]
  dim(plikes)
  1 - resample_misclassification(plikes, 1:K, m = k)
})
mean(aba_hats)
var(aba_hats)

## variance based on marginal stuff

Ksamp <- sample(bigK, K)
plikes <- plikes_pop[Ksamp, Ksamp]


bas <- resample_accs(plikes, k)
mean(bas)
var(bas)/K*k

## variance based on overlap covs

Ksamp <- sample(bigK, K)
plikes <- plikes_pop[Ksamp, Ksamp]
ol_covs <- sapply(1:k, function(i) cov_accs_overlap(plikes, k, i))
ol_covs

ms <- 1:k
sum(choose(K, ms) * choose(K - ms, k - ms) * choose(K - ms - k, k - ms) * ol_covs)/(choose(K, k)^2)




