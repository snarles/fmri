## develop a good basis to use
library(parallel)
source("extrapolation/ku_source.R")
source("extrapolation/gaussian_moments.R")
source("extrapolation/basis_source.R")

load("approximation/sim_large5_k5_raw.RData", verbose = TRUE)
mcc <- 20
sigma2_seq <- 0.005 * 1:50
mc.reps <- nrow(accsZ)
sigma2s <- rep(sigma2_seq, floor(mc.reps/length(sigma2_seq)))
facc <- t(sapply(1:50, function(ind) colMeans(accsZ[sigma2s==sigma2s[ind], ])))
facc_rep <- facc[match(sigma2s, sigma2_seq), ]
sqrt(colMeans((accsZ - facc_rep)^2))
kernel_sigma <- 0.5

(max.mu <- (qnorm(1- 1/(max(kref)^2))))

####
##  mass production
####
rmse_table <- function(accs_hat) {
  resids <- accs_hat - facc_rep
  resids2 <- resids^2
  biases <- t(sapply(1:50, function(ind) colMeans(resids[sigma2s==sigma2s[ind], ])))
  rmses <- t(sapply(1:50, function(ind) sqrt(colMeans(resids2[sigma2s==sigma2s[ind], ]))))
  list(biases = biases, rmses = rmses)
}

get_basis_mat2 <- function(max.mu, kernel_sigma, sub.rate = 1, mc.reps = 1e4) {
  (n_half <- ceiling(max.mu/kernel_sigma))
  (seq_half <- seq(0, max.mu, length.out = n_half))
  (kernel_mus <- c(rev(seq_half[-1]), seq_half, Inf))
  subseq <- seq(sub.rate, length(kref), by = sub.rate)
  Xmat <- sapply(kernel_mus, function(mu) gaussian_basis(mu, kernel_sigma^2, kref[subseq], mc.reps))
  Xtarg <- sapply(kernel_mus, function(mu) gaussian_basis(mu, kernel_sigma^2, Ktarg, mc.reps))
  list(Xmat = Xmat, Xtarg = Xtarg)
}

get_pred <- function(ind, Xmat, Xtarg, sub.rate = 1) {
  subseq <- seq(sub.rate, length(kref), by = sub.rate)
  nnls_fit <- nnls(Xmat, accs_subs[ind, subseq])
  as.numeric(Xtarg %*% nnls_fit$x)
}

sr <- 5
set1 <- get_basis_mat2(max.mu, 0.5, sub.rate = sr)
t1 <- proc.time()
preds <- mclapply(1:nrow(accsZ), get_pred, Xmat = set1$Xmat, Xtarg = set1$Xtarg, 
                  mc.cores = mcc, sub.rate = sr)
proc.time() - t1
accs_hat <- do.call(rbind, preds)

evalres <- rmse_table(accs_hat)
matplot(Ktarg, t(evalres$biases), type = "l")
matplot(Ktarg, t(evalres$rmses), type = "l")
apply(evalres$rmses, 2, max)

