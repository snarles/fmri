####
## BAYESIAN CLASSIFICATION
####

library(pracma)
library(magrittr)
library(MASS)
library(class)
source('transfer/source.R')

## PROBLEM PARAMETERS
n <- 30
p <- 20
sigma_x <- randn(2 * p, p) %>% {t(.) %*% .}
sigma_e <- randn(3 * p, p) %>% {t(.) %*% .}

## FORM MATRICES
blockmat0 <- rbind(
  cbind(sigma_x + sigma_e, sigma_x),
  cbind(sigma_x, sigma_x + sigma_e))
blockmat <- solve(blockmat0)
blockmat[1:p, 1:p] <- blockmat[1:p, 1:p] - solve(sigma_x + sigma_e)

## GENERATE DATA
mus <- mvrnorm(n, mu = rep(0, p), Sigma = sigma_x)
tr_set <- mus + mvrnorm(n, mu = rep(0, p), Sigma = sigma_e)
te_set <- mus + mvrnorm(n, mu = rep(0, p), Sigma = sigma_e)
te_ind <- sample(n, 1)

## NON BAYES APPROACH
a <- isqrtm(sigma_e)
te_cl <- knn(tr_set %*% a, te_set %*% a, 1:n)
(errs0 <- sum(te_cl != 1:n))

## APPLY BAYES
te_cl <- numeric(n)
for (te_ind in 1:n) {
  y <- te_set[te_ind, ]
  scores <- numeric(n)
  for (i in 1:n) {
    v <- c(tr_set[i, ], y)
    scores[i] <- t(v) %*% blockmat %*% v
  }
  te_cl[te_ind] <- which(scores == min(scores))
}
(errsB <- sum(te_cl != 1:n))
