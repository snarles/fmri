####
## COMPARE E BAYES WITH ML
###

library(magrittr)
library(pracma)
library(MASS)
library(glmnet)
source('transfer/source.R')
library(class)

## Generate parameters

n <- 100
p <- 10
q <- 20

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
Sigma_E_hat <- cov(resids)
Inv_Sigma_B <- diag(rep(1/lambdas, each = q))
EB_cov_B_vec <- solve(solve(Sigma_E_hat) %x% (t(X) %*% X) + Inv_Sigma_B)
f_EB_cov_mu <- function(xi, xj) {
  (diag(rep(1, p)) %x% t(xi)) %*% EB_cov_B_vec %*% (diag(rep(1, p)) %x% t(t(xj)))
}

# compare estimated sigma_B with true sigma_B
plot(s0s, lambdas)

## Generate new stimuli
L <- 100
n_te <- 20
X_te <- randn(L, q) 
i_chosen <- sample(L, n_te, TRUE)
y_star <- X_te[i_chosen, , drop = FALSE] %*% B0 +  mvrnorm(n_te, rep(0, p), Sigma_e)

## Maximum likelihood

mu_hat <- X_te %*% B_mu_EB
a <- isqrtm(Sigma_E_hat)
mle_cl <- knn(mu_hat %*% a, y_star %*% a, 1:L)
(mle_err <- sum(mle_cl != i_chosen))

## EM-bayes

mus <- list(L)
covs <- list(L)
cc <- rgb(0, 0, 1, alpha = 0.2)
pprobs <- matrix(0, n_te, L)
for (i in 1:L) {
  mus[[i]] <- X_te[i, , drop = FALSE] %*% B_mu_EB
  covs[[i]] <- f_EB_cov_mu(X_te[i, ], X_te[i, ])
  for (j in 1:n_te) {
    pprobs[j, i] <- -log(det(covs[[i]])) - 
      (y_star[j, , drop = FALSE] - mus[[i]]) %*% solve(covs[[i]]) %*% 
      t(y_star[j, , drop = FALSE] - mus[[i]])    
  }  
}
eb_cl <- apply(pprobs, 1, function(v) order(-v)[1])
(eb_err <- sum(eb_cl != i_chosen))
