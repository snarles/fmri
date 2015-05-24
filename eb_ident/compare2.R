####
## LARGE SCALE SIMULATION RESULTS
####

library(magrittr)
library(pracma)
library(MASS)
library(glmnet)
source('transfer/source.R')
source('eb_ident/eigenprism.R')
source('eb_ident/source.R')
library(class)

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
  ## Estimate parameters
  p_Bayes <- params_Bayes(X, Y, Sigma_e, Sigma_B, s0s)
  p_CV <- params_CV(X, Y, X, Y)
  p_EP <- params_EP(X, Y, X, Y)
  ## Generate new stimuli
  X_te <- randn(L, q) 
  i_chosen <- sample(L, n_te, TRUE)
  y_star <- X_te[i_chosen, , drop = FALSE] %*% B0 +  mvrnorm(n_te, rep(0, p), Sigma_e)
  ## error rates
  (bayes_err <- post_probs(X_te, y_star, i_chosen, p_Bayes)$err)
  (mle_err <- mle_cl(X_te, y_star, i_chosen, p_CV)$err)
  (mle_err_EP <- mle_cl(X_te, y_star, i_chosen, p_EP)$err)
  (eb_err <- post_probs(X_te, y_star, i_chosen, p_CV)$err)
  (ep_err <- post_probs(X_te, y_star, i_chosen, p_EP)$err)
  c(bayes_err, mle_err, eb_err, mle_err_EP, ep_err)/n_te
}


