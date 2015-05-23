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
B_mu <- matrix(0, q, p)
for (i in 1:p) {
  B_mu[, i] <- solve(t(X) %*% X + diag(rep(1/s0s[i], q))) %*% t(X) %*% Y[, i]    
}
Sigma_Bhat_vec <- solve(solve(Sigma_e) %x% (t(X) %*% X) + diag(1/diag(Sigma_B)))

## Estimate parameters

# using CV
B_mu_EB <- matrix(0, q, p) # placeholder for estimate
resids <- matrix(0, n, p)
lambdas <- numeric(p) # lambda for each column of B
for (i in 1:p) {
  res <- cv.glmnet(X, Y[, i], alpha = 0, intercept = FALSE, standardize = FALSE)
  pre <- predict(res, newx = X)
  resids[, i] <- pre - Y[, i]
  B_mu_EB[, i] <- as.numeric(coefficients(res))[-1]  
  lambdas[i] <- (sum(as.numeric(coefficients(res))[-1]^2)/q)  
}
lambdas <- pmax(lambdas, 1e-5)
Sigma_E_hat <- 0.5 * cov(resids) + 0.5 * diag(diag(cov(resids)))
Inv_Sigma_B <- diag(rep(1/lambdas, each = q))
EB_cov_B_vec <- solve(solve(Sigma_E_hat) %x% (t(X) %*% X) + Inv_Sigma_B)

# Using eigenprism
res_EP <- eigenprisms(X, Y)
filt_EP <- (res_EP$T2 > 0)
lambdas_EP <-  res_EP$T2[filt_EP]/q
p_EP <- sum(filt_EP)
B_mu_EP <- matrix(0, q, p_EP)
Y_EP <- Y[, filt_EP]
for (i in 1:p_EP) {
  B_mu_EP[, i] <- solve(t(X) %*% X + diag(rep(lambdas_EP[i], q))) %*% t(X) %*% Y_EP[, i]    
}
Yhat_EP <- X %*% B_mu_EP
Sigma_e_EP <- cov(Y_EP - Yhat_EP) %>% {0.5 * . + 0.5 * diag(diag(.))}
Sigma_B_EP <- diag(rep(1/lambdas_EP, each = q))
EP_cov_B_vec <- solve(solve(Sigma_e_EP) %x% (t(X) %*% X) + diag(1/diag(Sigma_B_EP)))

# compare estimated sigma_B with true sigma_B
#plot(s0s, lambdas)
plot(colSums(B0^2)/q, s0s)
points(colSums(B0^2)/q, lambdas, col = "blue")
points(colSums(B0^2)/q, res_EP$T2/q, col = "red")
abline(0, 1)

## Generate new stimuli
L <- 200
n_te <- 200
X_te <- randn(L, q) 
i_chosen <- sample(L, n_te, TRUE)
y_star <- X_te[i_chosen, , drop = FALSE] %*% B0 +  mvrnorm(n_te, rep(0, p), Sigma_e)
y_star_EP <- y_star[, filt_EP]

## Bayes
pp_bayes <- post_probs(Sigma_e, Sigma_Bhat_vec, B_mu, X_te, y_star)
bayes_cl <- apply(pp_bayes, 1, function(v) order(-v)[1])
(bayes_err <- sum(bayes_cl != i_chosen))

## Maximum likelihood
mu_hat <- X_te %*% B_mu_EB
a <- isqrtm(Sigma_E_hat)
mle_cl <- knn(mu_hat %*% a, y_star %*% a, 1:L)
(mle_err <- sum(mle_cl != i_chosen))

## Maximum likelihood with EP
mu_hat <- X_te %*% B_mu_EP
a <- isqrtm(Sigma_e_EP)
mle_cl_EP <- knn(mu_hat %*% a, y_star_EP %*% a, 1:L)
(mle_err_EP <- sum(mle_cl_EP != i_chosen))

## E-bayes
pp_eb <- post_probs(Sigma_E_hat, EB_cov_B_vec, B_mu_EB, X_te, y_star)
eb_cl <- apply(pp_eb, 1, function(v) order(-v)[1])
(eb_err <- sum(eb_cl != i_chosen))

## E-bayes with EP
pp_ep <- post_probs(Sigma_e_EP, EP_cov_B_vec, B_mu_EP, X_te, y_star_EP)
ep_cl <- apply(pp_ep, 1, function(v) order(-v)[1])
(ep_err <- sum(ep_cl != i_chosen))
