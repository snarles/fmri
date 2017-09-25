## develop a good basis to use
library(parallel)
source("extrapolation/ku_source.R")
source("extrapolation/gaussian_moments.R")

load("approximation/sim_large5_k5_raw.RData", verbose = TRUE)
mcc <- 20
sigma2_seq <- 0.005 * 1:50
mc.reps <- nrow(accsZ)
sigma2s <- rep(sigma2_seq, floor(mc.reps/length(sigma2_seq)))
facc <- t(sapply(1:50, function(ind) colMeans(accsZ[sigma2s==sigma2s[ind], ])))
facc_rep <- facc[match(sigma2s, sigma2_seq), ]
sqrt(colMeans((accsZ - facc_rep)^2))
#matplot(Ktarg, t(facc), type = "l", ylim = c(0, 1))

kernel_sigma <- 0.5

# plot(kref, accs_subs[5, ], type = "l")
# plot(Ktarg, accsZ[5, ], type = "l")

(max.mu <- (qnorm(1- 1/(max(kref)^2))))
(n_half <- ceiling(max.mu/kernel_sigma))
(seq_half <- seq(0, max.mu, length.out = n_half))
(kernel_mus <- c(rev(seq_half[-1]), seq_half, Inf))

Xmat <- sapply(kernel_mus, function(mu) gaussian_basis(mu, kernel_sigma^2, kref, 1e4))
Xtarg <- sapply(kernel_mus, function(mu) gaussian_basis(mu, kernel_sigma^2, Ktarg, 1e4))

#matplot(c(kref, Ktarg), rbind(Xmat, Xtarg), type = "l")


ind <- 20
true_acc <- colMeans(accsZ[sigma2s==sigma2s[ind], ])
accs_sub <- accs_subs[ind, ]

nnls_fit <- nnls(Xmat, accs_sub)
sum(nnls_fit$x)
acc_pred <- as.numeric(Xtarg %*% nnls_fit$x)

plot(c(kref, Ktarg), c(accs_sub, true_acc), type = "l", ylim = c(0, 1))
lines(c(kref, Ktarg), c(nnls_fit$fitted, acc_pred), col = "red")


####
##  mass production
####

get_basis_mat <- function(max.mu, kernel_sigma, mc.reps = 1e4) {
  (n_half <- ceiling(max.mu/kernel_sigma))
  (seq_half <- seq(0, max.mu, length.out = n_half))
  (kernel_mus <- c(rev(seq_half[-1]), seq_half, Inf))
  
  Xmat <- sapply(kernel_mus, function(mu) gaussian_basis(mu, kernel_sigma^2, kref, mc.reps))
  Xtarg <- sapply(kernel_mus, function(mu) gaussian_basis(mu, kernel_sigma^2, Ktarg, mc.reps))
  list(Xmat = Xmat, Xtarg = Xtarg)
}

get_pred <- function(ind, Xmat, Xtarg) {
  nnls_fit <- nnls(Xmat, accs_subs[ind, ])
  as.numeric(Xtarg %*% nnls_fit$x)
}

set1 <- get_basis_mat(max.mu, 0.25)
t1 <- proc.time()
preds <- mclapply(1:nrow(accsZ), get_pred, Xmat = set1$Xmat, Xtarg = set1$Xtarg, 
                  mc.cores = mcc)
proc.time() - t1


