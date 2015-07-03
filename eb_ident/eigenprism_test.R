## Regression with EigenPrism!
## Compare CV-ridge with eigenprism

library(magrittr)
library(pracma)
library(MASS)
library(glmnet)
source('eb_ident/eigenprism.R')

n <- 30
p <- 200
sigma_beta <- 10
sigma_y <- 20
sigma2_y <- sigma_y^2
W_X <- 2
Sigma_X <- 1/(W_X * p) * randn(W_X * p, p) %>% { t(.) %*% . }
X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma_X)
X_te <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma_X)
beta0 <- sigma_beta * rnorm(p)
T2_true <- Norm(beta0)^2
## True regression function
mu <- X %*% beta0
mu_te <- X_te %*% beta0

## Generate data
y <- mu + sigma_y * rnorm(n)

#rmse_y <- Norm(y - mu)

## Oracle Bayes
sigma2_b <- sigma_beta^2
lambda_bayes <- sigma2_y/sigma2_b
beta_bayes <- solve(t(X) %*% X + diag(rep(lambda_bayes, p)), t(X) %*% y)
yh_bayes <- X_te %*% beta_bayes
(rmse_bayes <- Norm(yh_bayes - mu_te))

## CV-glmnet
res <- cv.glmnet(X, y, intercept = FALSE, standardize = FALSE, alpha = 0)
beta_cv <- coef.cv.glmnet(res, s = "lambda.min", alpha = 0)[-1]
#yh_cv <- predict(res, X_te)
lambda_cv <- res$lambda.min * n /sd(y)
beta_cv2 <- solve(t(X) %*% X + diag(rep(lambda_cv, p)), t(X) %*% y)
f2(beta_cv2 - beta_cv)
yh_cv <- X_te %*% beta_cv
(rmse_cv <- Norm(yh_cv - mu_te))

## Eigenprism
ep <- eigenprism(X, y)
c(ep$T2, T2_true)
sigma2_b_hat <- ep$T2/p
c(sigma2_b_hat, sigma2_b)
eps <- eigenprism(X, y, sigma2 = TRUE)
(sigma2_y_hat <- eps$T2) ## noise estimation is too variable
sigma2_cv <- with(res, cvm[which(lambda == res$lambda.1se)])/n
(lambda_ep <- sigma2_cv/sigma2_b_hat)
beta_ep <- solve(t(X) %*% X + diag(rep(lambda_ep, p)), t(X) %*% y)
yh_ep <- X_te %*% beta_ep
(rmse_ep <- Norm(yh_ep - mu_te))

"Lambdas"
c(Bayes = lambda_bayes, CV = lambda_cv, EP = lambda_ep)
"Norm(pred - truth)"
c(Bayes = rmse_bayes, CV = rmse_cv, EP = rmse_ep)
