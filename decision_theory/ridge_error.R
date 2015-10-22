library(MASS)
library(pracma)
library(parallel)

f2 <- function(x, y=0) sum((x-y)^2)

ridge_error_trial <- function(bt, X, lambda) {
  #p <- dim(Sigma)[1]
  #X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  y <- X %*% bt + rnorm(dim(X)[1])
  btr <- solve(t(X) %*% X + lambda * eye(dim(X)[2]), t(X) %*% y)
  f2(bt, btr)
}

ridge_error_trials <- function(bt, X, lambda,
                               mc.reps = 100, mc.cores = 7) {
  res <- mclapply(1:mc.reps, function(i) {
    set.seed(i)
    ridge_error_trial(bt, X, lambda)
  }, mc.cores = mc.cores)
  mean(unlist(res))
}


ridge_error_theory <- function(bt, X, lambda) {
  res <- svd(X); V <- res$v; d <- res$d; p <- dim(X)[2]
  vbt <- t(V) %*% bt
  diags <- (lambda/(d^2 + lambda)) - 1  
  diags2 <- (lambda/(d^2 + lambda))^2 - 1
  bias_term <- f2(bt) + sum(diags2 * vbt^2)
  var_term <- 1/lambda * (p + sum(diags)) - 1/lambda * (p + sum(diags2))
  temp <- solve(t(X) %*% X + lambda * eye(p), t(X))
  var_term <- sum(diag(temp %*% t(temp)))
  bias_term + var_term
}

ridge_error_theory_ <- function(bt, X) {
  res <- svd(X); V <- res$v; d <- res$d; p <- dim(X)[2]
  vbt <- t(V) %*% bt
  fbt <- f2(bt)
  ff <- function(lambda) {
    diags <- (lambda/(d^2 + lambda)) - 1  
    diags2 <- (lambda/(d^2 + lambda))^2 - 1
    f2(bt) + sum(diags2 * vbt^2) + 
      1/lambda * (p + sum(diags)) - 1/lambda * (p + sum(diags2))
  }
  ff
}

ridge_oracle_error <- function(bt, X, bounds = c(1e-5, 1e5)) {
  ff <- ridge_error_theory_(bt, X)
  res <- optimise(ff, bounds)
  res
}

# bt <- rnorm(10); X <-randn(5, 10); lambda <- 1 
# ridge_error_trials(bt, X, lambda, mc.reps = 2000)
# ridge_error_theory(bt, X, lambda)
# sapply(seq(0.1, 1, 0.1), ridge_error_theory_(bt, X))
# optimise(ridge_error_theory_(bt, X), c(1e-5, 1e5))