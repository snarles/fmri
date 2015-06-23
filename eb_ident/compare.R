library(magrittr)
library(pracma)
library(MASS)
library(glmnet)
source('transfer/source.R')
source('eb_ident/eigenprism.R')
source('eb_ident/source.R')
library(class)

#set.seed(1)
pars <- list(n=30, p= 60, q= 60, W_X= 2, s_e= 10, W_e= 2, L= 100, n_te= 10)
zattach(pars)


Sigma_X <- 1/(W_X * q) * randn(W_X * q, q) %>% { t(.) %*% . }
Sigma_e <- s_e/(W_e * p) * randn(W_e * p, p) %>% { t(.) %*% . }
s0s <- 0* rnorm(p)^2 + 1
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
p_EPa <- params_EP(X, Y, X, Y, TRUE)
p_sd <- p_EP_2_sd(X, Y, p_EPa)

p_Bayes_wrong_e <- params_Bayes(X, Y, p_CV$Sigma_e, Sigma_B, s0s)
#p_Bayes_wrong_e <- p_Bayes
#p_Bayes_wrong_e$Sigma_e <- p_CV$Sigma_e

## Generate new stimuli
X_te <- randn(L, q)
i_chosen <- sample(L, n_te, TRUE)
y_star <- X_te[i_chosen, , drop = FALSE] %*% B0 +  mvrnorm(n_te, rep(0, p), Sigma_e)
## error rates
(bayes_err <- post_probs(X_te, y_star, i_chosen, p_Bayes)$err)

(bayes_wrong_e_err <- post_probs(X_te, y_star, i_chosen, p_Bayes_wrong_e)$err)

(mle_err <- mle_cl(X_te, y_star, i_chosen, p_CV)$err)
(mle_err_EP <- mle_cl(X_te, y_star, i_chosen, p_EP)$err)
(mle_oracle_err <- mle_cl(X_te, y_star, i_chosen, p_Bayes)$err)

(eb_err <- post_probs(X_te, y_star, i_chosen, p_CV)$err)
(ep_err <- post_probs(X_te, y_star, i_chosen, p_EP)$err)
(epa_err <- post_probs(X_te, y_star, i_chosen, p_EPa)$err)
(sd_err <- post_probs(X_te, y_star, i_chosen, p_sd)$err)
