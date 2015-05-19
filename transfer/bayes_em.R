####
## EMPIRICAL BAYESIAN CLASSIFICATION
####

library(pracma)
library(magrittr)
library(MASS)
library(class)
source('transfer/source.R')

## PROBLEM PARAMETERS
n <- 100
p <- 50
sigma_mu <- randn(2 * p, p) %>% {t(.) %*% .}
sigma_e <- randn(3 * p, p) %>% {t(.) %*% .}
sigma_y <- sigma_mu + sigma_e
mu_mu <- rnorm(p)
sigma_compound <- rbind(cbind(sigma_y, sigma_mu), cbind(sigma_mu, sigma_y))

## BAYES RULE
bayes_mat <- solve(sigma_compound)
bayes_mat[1:p, 1:p] <- bayes_mat[1:p, 1:p] - solve(sigma_y)

## GENERATE DATA
mus <- mvrnorm(n, mu = mu_mu, Sigma = sigma_mu)
tr_set <- mus + mvrnorm(n, mu = rep(0, p), Sigma = sigma_e)
te_set1 <- mus + mvrnorm(n, mu = rep(0, p), Sigma = sigma_e)
te_set2 <- mus + mvrnorm(n, mu = rep(0, p), Sigma = sigma_e)

## ESTIMATE COVARIANCES

te_set <- te_set1
est_parameters <- function(tr_set, te_set, set_in_env = FALSE) {
  n <- dim(te_set)[1]
  mu_est <- colMeans(rbind(tr_set, te_set))
  cent_set <- t(t(rbind(tr_set, te_set)) - mu_est)
  cov_y_est <- t(cent_set) %*% cent_set / (2 * n)
  cent_tr <- t(t(tr_set) - mu_est)
  cent_te <- t(t(te_set) - mu_est)
  cov_mu_est <- (t(cent_tr) %*% cent_te + t(cent_te) %*% cent_tr) / (2 * n)
  if (set_in_env) {
    mu_est <<- mu_est
    cov_y_est <<- cov_y_est
    cov_mu_est <<- cov_mu_est
    return()
  }
  list(mu_est = mu_est, cov_y_est = cov_y_est, cov_mu_est = cov_mu_est)
}

est_parameters(tr_set, te_set1, TRUE)

plot(mu_est, mu_mu)
plot(cov_y_est, sigma_y); abline(0, 1)
plot(cov_mu_est, sigma_mu); abline(0, 1)
svd(cov_y_est)$d
svd(cov_mu_est)$d
cov_compound <- rbind(cbind(cov_y_est, cov_mu_est), cbind(cov_mu_est, cov_y_est))
svd(cov_compound)$d

## EMPRIRICAL BAYES RULE

Ebayes_mat <- solve(cov_compound)
Ebayes_mat[1:p, 1:p] <- Ebayes_mat[1:p, 1:p] - solve(cov_y_est)

####
## PERFORMANCE OF EMPIRICAL BAYES VS BAYES
####

classify <- function(tr_set, te_set, mat, mu_est) {
  te_cl <- numeric(n)
  for (te_ind in 1:n) {
    y <- te_set[te_ind, ]
    scores <- numeric(n)
    for (i in 1:n) {
      v <- c(tr_set[i, ], y) - mu_est
      scores[i] <- t(v) %*% mat %*% v
    }
    te_cl[te_ind] <- which(scores == min(scores))
  }
  te_cl
}

## BAYES
bayes_cl <- classify(tr_set, te_set2, bayes_mat, mu_mu)
sum(bayes_cl == 1:n)

## EM-BAYES
Ebayes_cl <- classify(tr_set, te_set2, Ebayes_mat, mu_est)
sum(Ebayes_cl == 1:n)
