gcv_objective_ <- function(X, y) {
  res <- svd(X)
  U <- res$u
  d <- res$d
  xi <- as.numeric(t(U) %*% y)  
  yn <- sum(y^2)
  n <- dim(y)[1]
  ff <- function(lambda) {
    (yn + sum(xi^2 * ((d^2/(d^2 + lambda))^2 - 2 * (d^2/(d^2 + lambda)))))/
      ( n - sum(d^2/(d^2 + lambda)))
  }  
  ff
}

# library(pracma)
# 
# f2 <- function(x, y=0) sum((x-y)^2)
# 
# n <- 30; p <- 200;
# X <- randn(n, p)
# bt <- rnorm(p)
# y <- as.numeric(X %*% bt + rnorm(n))
# res <- svd(X)
# U <- res$u
# d <- res$d
# xi <- as.numeric(t(U) %*% y)
# 
# ## verify some formulas
# lambda <- 2
# yhat <- X %*% solve(t(X) %*% X + lambda * eye(p), t(X) %*% y)
# yhat2 <- U %*% diag(d^2/(d^2 + lambda)) %*% t(U) %*% y
# (trr <- sum(diag(X %*% solve(t(X) %*% X + lambda * eye(p), t(X)))))
# sum(d^2/(d^2 + lambda))
# f2(y - yhat)
# f2(y) + sum(xi^2 * ((d^2/(d^2 + lambda))^2 - 2 * (d^2/(d^2 + lambda))))
# 
# f2(yhat, yhat2)


