library(magrittr)
library(MASS)
source("decision_theory//gcv.R")
source("decision_theory//ridge_error.R")

risk_ratio_exp <- function(p, n, alpha, vs_X, ws_X,
                           mc.reps = 1e4, mc.cores = 7) {
  U <- svd(randn(p, p))$u
  d_X <- sample(vs_X, p, prob=ws_X, replace = TRUE)
  Sigma_X <- t(U) %*% diag(d_X) %*% U
  X <- mvrnorm(n, rep(0, p), Sigma_X)
  bt <- rnorm(p) %>% {sqrt(alpha) * ./sqrt(f2(.))}
  
  ridge_error_theory(bt, X, 0)
  oracle_res <- ridge_oracle_error(bt, X)
  oracle_lambda <- oracle_res[[1]]
  oracle_error <- oracle_res[[2]]
  
  gcv_res <- gcv_trials(bt, X, obj = gcv_objective_,
                        mc.reps = mc.reps, mc.cores = mc.cores)
  gcv_lambda <- median(gcv_res[, 1])
  gcv_error <- mean(gcv_res[, 2])
  press_res <- gcv_trials(bt, X, obj = press_objective_,
                        mc.reps = mc.reps, mc.cores = mc.cores)
  press_lambda <- median(press_res[, 1])
  press_error <- mean(press_res[, 2])  
  
  print(paste("lambdas: oracle", oracle_lambda, "gcv", gcv_lambda,
              "press", press_lambda))
  print(paste("error: oracle", oracle_error, "gcv", gcv_error, 
              "ratio", gcv_error/oracle_error, "press", press_error,
              "ratio", press_error/oracle_error))
  list(X = X, bt = bt, Sigma_X = Sigma_X, oracle_lambda = oracle_lambda,
       oracle_error = oracle_error, gcv_res = gcv_res, gcv_lambda = gcv_lambda,
       gcv_error = gcv_error, press_res = press_res, press_lambda = press_lambda,
       press_error = press_error,
       p = p, n = n, alpha = alpha, vs_X = vs_X, ws_X = ws_X)
}

res <- risk_ratio_exp(
  p = 10,
  n = 5,
  alpha = 0,
  vs_X = c(1.0, 0.0),
  ws_X = c(1.0, 0.0),
  mc.reps = 1e1
)

lineId::zattach(res)
ff <- ridge_error_theory_(bt, X)
oracle_lambda
ff(oracle_lambda)
ridge_error_trials(bt, X, oracle_lambda) %>% summary

y <- X %*% bt + rnorm(n)
ff <- press_objective_(X, y)
vals <- 10^seq(-5, 5, by = .1)
plot(log(vals), sapply(vals, ff), type = 'l')
