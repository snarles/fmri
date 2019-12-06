library(pracma)
library(MASS)
library(class)
library(parallel)
library(FNN)
library(lineId) ## use devtools::install('lineId')
#source('extrapolation/kay_method.R')
#source('fingerprint/mi_est_pipeline.R')


mi_simulation <- function(K, d, p, sigma2) {
  xs <- randn(K, d)
  amat <- randn(d, p)
  ys <- (xs + sqrt(sigma2) * randn(K, d)) %*% amat
  ys2 <- (xs + sqrt(sigma2) * randn(K, d)) %*% amat
  
  # true MI between X and Y, I(X; Y)
  rho_true <- 1/sqrt(1 + sigma2)
  mi_true <- d * (-.5 * log(1 - rho_true^2))
  mi_true
  
  # scores matrix
  pmat <- -pdist2(ys, ys2)
  resample_misclassification(pmat, 1:K)
  
  K <- dim(pmat)[1]
  empirical_acc <- 1- resample_misclassification(pmat, 1:K)
  
  mi_grid <- 3 * (1:600)/100
  pi_ests <- piK(sqrt(2 * mi_grid), K)
  mi_xx_est <- max(mi_grid[(1-pi_ests) < empirical_acc ])
  mi_xx_est
  
  component_mis <- sapply(1:p, function(i) mutinfo(ys[, i], ys2[, i]))
  sum(component_mis)
  
  ratio <- sum(component_mis)/mi_xx_est
  est_d <- p/ratio
  est_d
  
  # mi per dimension
  mi_per_dim <- mi_xx_est/est_d
  # invert to rho since -1/2 * log(1 - rho^4) = I(X;X')
  rho_est <- (1 - exp(-2 * mi_per_dim)) ^ (.25)
  mi_xy_est <- -est_d/2 * log(1 - rho_est^2)
  c(xx=mi_xx_est, est=mi_xy_est, true=mi_true, d_est = est_d)
}


nreps <- 40
dpars <- c(5, 10, 15, 20, 25, 30, 35, 40, 45)
exp_params <- cbind(80, rep(dpars, each = nreps), 100, 0.5)
results <- apply(exp_params, 1, function(x) mi_simulation(x[1], x[2], x[3], x[4]))
saveRDS(results, 'fingerprint/mi_xx_to_xy_sim.rds')

#results
results2 <- array(results, dim=c(4, nreps, length(dpars)))
meds <- apply(results2, c(1, 3), median)
plot(results[3, ], results[2, ], ylim = c(0, 30), xlab='true I(X; Y)', ylab='est', col = 'red')
points(results[3, ], results[1, ], col='blue')
lines(meds[3, ], meds[2, ], col = 'red')
lines(meds[3, ], meds[1, ], col = 'blue')
abline(0, 1)
legend(3, 20, c(expression(hat(I)(X,Y)), "median", expression(hat(I)(X, X^I)), "median"), 
       col=c('red', 'red', 'blue', 'blue'), pch=c('o', NA, 'o', NA), lty=c(NA, 1, NA, 1))

