## Regression with EigenPrism!
## Compare CV-ridge with eigenprism

library(magrittr)
library(pracma)
library(MASS)
library(glmnet)
source('transfer/source.R')
source('eb_ident/eigenprism.R')

n <- 30
p <- 200
sigma_beta <- 2
sigma_y <- 10
W_X <- 2
Sigma_X <- 1/(W_X * p) * randn(W_X * p, p) %>% { t(.) %*% . }
X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma_X)
X <- scale(X, center = TRUE, scale = TRUE)
beta0 <- sigma_beta * rnorm(p)
## True regression function
mu <- X %*% beta0

## Generate data
y <- mu + sigma_y * rnorm(n)

## CV-glmnet
res <- cv.glmnet(X, y, intercept = FALSE)
beta_cv <- coef.cv.glmnet(res, s = "lambda.1se")
yh_cv <- predict(res, X)
Norm(yh_cv - mu)

## Eigenprism
ep <- eigenprism(X, y)
ep$T2
beta_ep <- solve(t(X) %*% X + diag(rep(1/ep$T2, p)), t(X) %*% y)
yh_ep <- X %*% beta_ep
Norm(yh_ep - mu)
