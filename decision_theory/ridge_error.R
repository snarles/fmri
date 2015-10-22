library(MASS)
library(pracma)
library(parallel)

f2 <- function(x, y=0) sum((x-y)^2)

ridge_error_trial <- function(bt, X, lambda) {
  #p <- dim(Sigma)[1]
  #X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  y <- X %*% bt + rnorm(dim(X)[1])
  if (lambda == 0) {
    btr <- ginv(X) %*% y
  } else {
    btr <- solve(t(X) %*% X + lambda * eye(dim(X)[2]), t(X) %*% y)    
  }
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
  diags2 <- (lambda/(d^2 + lambda))^2 - 1
  bias_term <- f2(bt) + sum(diags2 * vbt^2)
  var_term <- sum(d^2/(d^2 + lambda)^2) 
  bias_term + var_term
}

ridge_error_theory_ <- function(bt, X) {
  res <- svd(X); V <- res$v; d <- res$d; p <- dim(X)[2]
  vbt <- t(V) %*% bt
  fbt <- f2(bt)
  ff <- function(lambda) {
    diags2 <- (lambda/(d^2 + lambda))^2 - 1
    f2(bt) + sum(diags2 * vbt^2) + 
      sum(d^2/(d^2 + lambda)^2)  
  }
  ff
}

ridge_oracle_error <- function(bt, X, bounds = c(1e-5, 1e5)) {
  ff <- ridge_error_theory_(bt, X)
  res <- optimise(ff, bounds)
  res
}

ridge_error_random_I <- function(alpha2, gamma, lambda, n = 1e3, naive = FALSE) {
  p <- floor(gamma * n)
  bt <- rnorm(p)
  bt <- sqrt(alpha2) * bt/sqrt(sum(bt^2))
  X <- randn(n, p)
  if (naive) {
    y <- X %*% bt + rnorm(n)
    bth <- solve(t(X) %*% X + lambda * eye(p), t(X) %*% y)
    x_star <- rnorm(p)
    y_star <- sum(x_star * bt) + rnorm(1)
    return((y_star - sum(x_star * bth))^2)
  } else {
    return(1 + ridge_error_theory(bt, X, lambda))
  }
}

# bt <- rnorm(10); X <-randn(5, 10); lambda <- 1 
# ridge_error_trials(bt, X, lambda, mc.reps = 2000)
# ridge_error_theory(bt, X, lambda)
# sapply(seq(0.1, 1, 0.1), ridge_error_theory_(bt, X))
# optimise(ridge_error_theory_(bt, X), c(1e-5, 1e5))