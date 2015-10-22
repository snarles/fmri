library(MASS)
library(pracma)
library(parallel)

f2 <- function(x, y=0) sum((x-y)^2)

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

gcv_trial <- function(bt, X, mu = X %*% bt, bounds = c(1e-5, 1e5),
                      obj = gcv_objective_) {
  y <- mu + rnorm(dim(X)[1])
  ff <- obj(X, y)
  res <- optimise(ff, bounds)
  lambda <- res$minimum
  bt_g <- solve(t(X) %*% X + lambda * eye(dim(X)[2]), t(X) %*% y)
  c(lambda = lambda, error = f2(bt, bt_g))
}

gcv_trials <- function(bt, X, bounds = c(1e-5, 1e5),
                       obj = gcv_objective_,
                       mc.reps = 100, mc.cores = 7) {
  mu <- X %*% bt
  res <- mclapply(1:mc.reps, function(i) {
    set.seed(i)
    gcv_trial(bt, X, mu, bounds, obj)
  }, mc.cores = mc.cores)
  do.call(rbind, res)
}


press_objective_ <- function(X, y) {
  res <- svd(X)
  U <- res$u
  d <- res$d
  #xi <- as.numeric(t(U) %*% y)  
  #yn <- sum(y^2)
  n <- dim(y)[1]
  ff <- function(lambda) {
    #amat <- X %*% solve(t(X) %*% X + lambda * eye(dim(X)[2]), t(X))
    amat <- U %*% diag(d^2/(d^2 + lambda)) %*% t(U)
    aks <- diag(amat)
    resids <- y - amat %*% y
    sum(resids^2 * (1/(1-aks))^2)
  }  
  ff
}

press_objective_naive <- function(X, y, lambda) {
  n <- dim(X)[1]; p <- dim(X)[2]
  errs <- sapply(1:n, function(i) {
    Xtr <- X[-i, ]; ytr <- y[-i]
    bth <- solve(t(Xtr) %*% Xtr + lambda * eye(p), t(Xtr) %*% ytr)
    yh <- X[i, , drop = FALSE] %*% bth 
    (yh - y[i])^2
  })
  sum(errs)
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
# invmat <- solve(t(X) %*% X + lambda * eye(p))
# invmat2 <- 1/(lambda) * (eye(p) + V %*% diag((lambda/(d^2 + lambda)) - 1) %*% t(V))
# f2(invmat, invmat2)
# sqmat <- solve(t(X) %*% X + lambda * eye(p)) %*% 
#   solve(t(X) %*% X + lambda * eye(p))
# sqmat2 <- 1/(lambda^2) * (eye(p) + V %*% diag((lambda/(d^2 + lambda))^2 - 1) %*% t(V))
# f2(sqmat, sqmat2)
# 
# temp <- solve(t(X) %*% X + lambda * eye(p), t(X))
# varmat <- temp %*% t(temp)
# varmat2 <- invmat - lambda * sqmat
# f2(varmat, varmat2)
