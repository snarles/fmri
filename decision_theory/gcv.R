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
# err <- rnorm(n)
# y <- as.numeric(X %*% bt + err)
# res <- svd(X)
# U <- res$u; V <- res$v; d <- res$d
# dim(V)
# xi <- as.numeric(t(U) %*% y)
# 
# ## verify some formulas
# lambda <- 2
# yhat <- X %*% solve(t(X) %*% X + lambda * eye(p), t(X) %*% y)
# yhat2 <- U %*% diag(d^2/(d^2 + lambda)) %*% t(U) %*% y
# bhat <-solve(t(X) %*% X + lambda * eye(p), t(X) %*% y)
# bhat2 <- bt - lambda * solve(t(X) %*% X + lambda * eye(p), bt) + 
#   solve(t(X) %*% X + lambda * eye(p), t(X) %*% err)
# bias <- bt - solve(t(X) %*% X + lambda * eye(p), t(X) %*% X %*% bt)
# bias2 <- (V %*% diag(lambda/(d^2 + lambda) - 1) %*% t(V) + eye(p)) %*% bt
# f2(bias, bias2)
# f2(bhat, bhat2)
# (trr <- sum(diag(X %*% solve(t(X) %*% X + lambda * eye(p), t(X)))))
# sum(d^2/(d^2 + lambda))
# f2(y - yhat)
# f2(y) + sum(xi^2 * ((d^2/(d^2 + lambda))^2 - 2 * (d^2/(d^2 + lambda))))
# 
# f2(yhat, yhat2)
# sqmat <- solve(t(X) %*% X + lambda * eye(p)) %*% 
#   solve(t(X) %*% X + lambda * eye(p))
# sqmat2 <- 1/(lambda^2) * (eye(p) + V %*% diag((lambda/(d^2 + lambda))^2 - 1) %*% t(V))
# f2(sqmat, sqmat2)


