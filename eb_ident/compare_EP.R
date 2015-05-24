####
## COMPARE E BAYES (+ with eigenprism) WITH ML
###

library(magrittr)
library(pracma)
library(MASS)
library(glmnet)
library(class)
source('transfer/source.R')
source('eb_ident/source.R')
source('eb_ident/eigenprism.R')

## Generate parameters

n <- 30
p <- 60
q <- 60

Sigma_X <- 1/(2 * q) * randn(2 * q, q) %>% { t(.) %*% . }
Sigma_e <- 10/(2 * p) * randn(2 * p, p) %>% { t(.) %*% . }
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
ps_bagged <- mle_bagged(X, Y)

# compare estimated sigma_B with true sigma_B
#plot(s0s, lambdas)
plot(colSums(B0^2)/q, s0s)
points(colSums(B0^2)/q, p_CV$lambdas, col = "blue")
points(colSums(B0^2)/q, p_EP$res_EP$T2/q, col = "red")
abline(0, 1)

## Generate new stimuli
L <- 200
n_te <- 200
X_te <- randn(L, q) 
i_chosen <- sample(L, n_te, TRUE)
y_star <- X_te[i_chosen, , drop = FALSE] %*% B0 +  mvrnorm(n_te, rep(0, p), Sigma_e)

## error rates
(bayes_err <- post_probs(X_te, y_star, i_chosen, p_Bayes)$err)
(mle_err <- mle_cl(X_te, y_star, i_chosen, p_CV)$err)
(mle_err_EP <- mle_cl(X_te, y_star, i_chosen, p_EP)$err)
(eb_err <- post_probs(X_te, y_star, i_chosen, p_CV)$err)
(ep_err <- post_probs(X_te, y_star, i_chosen, p_EP)$err)
(bagged_err <- bagged_cl(X_te, y_star, i_chosen, ps_bagged)$err)
