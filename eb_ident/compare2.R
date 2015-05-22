####
## LARGE SCALE SIMULATION RESULTS
####

library(magrittr)
library(pracma)
library(MASS)
library(glmnet)
source('transfer/source.R')
library(class)

post_probs <- function(Sigma_e, Sigma_vec_B, B_hat, X, Y) {
  L <- dim(X)[1]
  n <- dim(Y)[1]
  mus <- list(L)
  covs <- list(L)
  pprobs <- matrix(0, n, L)
  p <- dim(Y)[2]
  for (i in 1:L) {
    mus[[i]] <- X[i, , drop = FALSE] %*% B_hat
    covs[[i]] <- (diag(rep(1, p)) %x% t(X[i, ])) %*% 
      Sigma_vec_B %*% (diag(rep(1, p)) %x% t(t(X[i, ]))) + Sigma_e
    for (j in 1:n) {
      pprobs[j, i] <- -log(det(covs[[i]])) - 
        (Y[j, , drop = FALSE] - mus[[i]]) %*% solve(covs[[i]]) %*% 
        t(Y[j, , drop = FALSE] - mus[[i]])
    }  
  }
  pprobs
}

#' Parameters
#' @param n Training set size
#' @param p Dimension of response
#' @param q Dimension of feature
#' @param W_X Wishart parameter for feature (determines collinearity)
#' @param s_e Scale of noise
#' @param W_e Wishart parameter for noise
#' @param L Number of test classes
#' @param n_te Size of test set
run_exp <- function(n, p, q, W_X, s_e, W_e, 
                    L, n_te) {
  ## Generate parameters
  Sigma_X <- 1/(W_X * q) * randn(W_X * q, q) %>% { t(.) %*% . }
  Sigma_e <- s_e/(W_e * p) * randn(W_e * p, p) %>% { t(.) %*% . }
  s0s <- rnorm(p)^2
  Sigma_B <- diag(rep(s0s, each = q))
  Bvec <- sqrt(Sigma_B) %*% rnorm(p * q)
  B0 <- matrix(Bvec, q, p)
  ## Generate data
  X <- mvrnorm(n, mu = rep(0, q), Sigma = Sigma_X)
  Y <- X %*% B0 + mvrnorm(n, mu = rep(0, p), Sigma = Sigma_e)
  Sigma_Bhat_vec <- solve(solve(Sigma_e) %x% (t(X) %*% X) + diag(1/diag(Sigma_B)))
  B_mu <- matrix(0, q, p)
  for (i in 1:p) {
    B_mu[, i] <- solve(t(X) %*% X + diag(rep(1/s0s[i], q))) %*% t(X) %*% Y[, i]    
  }
  ## Estimate parameters
  B_mu_EB <- matrix(0, q, p) # placeholder for estimate
  resids <- matrix(0, n, p)
  lambdas <- numeric(p) # lambda for each column of B
  for (i in 1:p) {
    res <- cv.glmnet(X, Y[, i], alpha = 0, intercept = FALSE, 
                     standardize = FALSE, nfolds = 5)
    pre <- predict(res, newx = X)
    resids[, i] <- pre - Y[, i]
    B_mu_EB[, i] <- as.numeric(coefficients(res))[-1]  
    lambdas[i] <- (sum(as.numeric(coefficients(res))[-1]^2)/q)  
  }
  lambdas[lambdas < 1e-3] <- 1e-3
  Sigma_E_hat <- 0.5 * cov(resids) + 0.5 * diag(diag(cov(resids)))
  Inv_Sigma_B <- diag(rep(1/lambdas, each = q))
  EB_cov_B_vec <- solve(solve(Sigma_E_hat) %x% (t(X) %*% X) + Inv_Sigma_B)
  ## Generate new stimuli
  X_te <- randn(L, q) 
  i_chosen <- sample(L, n_te, TRUE)
  y_star <- X_te[i_chosen, , drop = FALSE] %*% B0 +  mvrnorm(n_te, rep(0, p), Sigma_e)
  ## Bayes
  pp_bayes <- post_probs(Sigma_e, Sigma_Bhat_vec, B_mu, X_te, y_star)
  bayes_cl <- apply(pp_bayes, 1, function(v) order(-v)[1])
  bayes_err <- sum(bayes_cl != i_chosen)/n_te
  ## Maximum likelihood
  mu_hat <- X_te %*% B_mu_EB
  a <- isqrtm(Sigma_E_hat)
  mle_cl <- knn(mu_hat %*% a, y_star %*% a, 1:L)
  mle_err <- sum(mle_cl != i_chosen)/n_te
  ## E-bayes
  pp_eb <- post_probs(Sigma_E_hat, EB_cov_B_vec, B_mu_EB, X_te, y_star)
  eb_cl <- apply(pp_eb, 1, function(v) order(-v)[1])
  eb_err <- sum(eb_cl != i_chosen)/n_te
  c(bayes_err, mle_err, eb_err)
}


