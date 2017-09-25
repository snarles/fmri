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
rmse_table <- function(accs_hat) {
  resids <- accs_hat - facc_rep
  resids2 <- resids^2
  biases <- t(sapply(1:50, function(ind) colMeans(resids[sigma2s==sigma2s[ind], ])))
  rmses <- t(sapply(1:50, function(ind) sqrt(colMeans(resids2[sigma2s==sigma2s[ind], ]))))
  list(biases = biases, rmses = rmses)
}

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

set1 <- get_basis_mat(max.mu, 0.5)
t1 <- proc.time()
preds <- mclapply(1:nrow(accsZ), get_pred, Xmat = set1$Xmat, Xtarg = set1$Xtarg, 
                  mc.cores = mcc)
proc.time() - t1
accs_hat <- do.call(rbind, preds)

evalres <- rmse_table(accs_hat)
matplot(Ktarg, t(evalres$biases), type = "l")
matplot(Ktarg, t(evalres$rmses), type = "l")
apply(evalres$rmses, 2, max)

## cross-validating the bandwidth??

bdwids <- seq(0.1, 1, by = 0.1)
basis_sets <- lapply(bdwids, function(x) get_basis_mat(max.mu, x))

get_fit <- function(ind, Xmat) {
  nnls_fit <- nnls(Xmat, accs_subs[ind, ])
  sqrt(mean(nnls_fit$residuals^2))
}

bdwid_cv_curve <- function(accs_sub, basis_sets) {
  ans <- numeric()
  for (j in 1:length(basis_sets)) {
    set1 <- basis_sets[[j]]
    nnls_fit <- nnls(Xmat[-length(accs_sub), ], accs_sub[-length(accs_sub)])
    pred <- sum(Xmat[length(accs_sub), ] * nnls_fit$x)
    ans[j] <- pred - accs_sub[length(accs_sub)]
  }
  ans
}

bdwid_cv_curve <- function(accs_sub, basis_sets, cv.frac = 0.5) {
  ans <- numeric()
  ntr <- length(accs_sub)
  tr.inds <- 1:floor(ntr * cv.frac)
  for (j in 1:length(basis_sets)) {
    set1 <- basis_sets[[j]]
    Xmat <- set1$Xmat
    nnls_fit <- nnls(Xmat[tr.inds, ], accs_sub[tr.inds])
    pred <- sum(Xmat[-tr.inds, ] * nnls_fit$x)
    ans[j] <- sqrt(mean((pred - accs_sub[-tr.inds])^2))
  }
  ans
}

bdwid_fit_curve <- function(accs_sub, basis_sets) {
  ans <- numeric()
  ntr <- length(accs_sub)
  for (j in 1:length(basis_sets)) {
    set1 <- basis_sets[[j]]
    Xmat <- set1$Xmat
    nnls_fit <- nnls(Xmat, accs_sub)
    ans[j] <- sqrt(mean(nnls_fit$residuals^2))
  }
  ans
}


fsizes <- sapply(basis_sets, function(set1) {
  mean(unlist(mclapply(1:nrow(accsZ), get_fit, Xmat = set1$Xmat, 
           mc.cores = mcc))^2)
})
plot(bdwids, fsizes)

rank(abs(bdwid_cv_curve(accs_subs[157, ], basis_sets)))
plot(abs(bdwid_cv_curve(accs_subs[157, ], basis_sets)))
plot(bdwid_cv_curve(accs_subs[22, ], basis_sets, cv.frac = 0.9))


plot(bdwid_fit_curve(accs_subs[22, ], basis_sets))
plot(bdwid_fit_curve(accs_subs[28, ], basis_sets)); abline(h = 1/max(kref))

####
##  method for choosing bandwidth and predicting
####

bdwid_select_predict <- function(accs_sub, basis_sets, kref) {
  curve <- bdwid_fit_curve(accs_sub, basis_sets)
  basis_ind <- max(which(curve < min(curve) + 1/max(kref)))
  nnls_fit <- nnls(basis_sets[[basis_ind]]$Xmat, accs_subs[ind, ])
  as.numeric(basis_sets[[basis_ind]]$Xtarg %*% nnls_fit$x)
}

t1 <- proc.time()
preds <- mclapply(1:nrow(accsZ), function(ind) bdwid_select_predict(accs_subs[ind, ], basis_sets, kref),
                  mc.cores = mcc)
proc.time() - t1
accs_hat <- do.call(rbind, preds)

evalres <- rmse_table(accs_hat)
matplot(Ktarg, t(evalres$biases), type = "l")
matplot(Ktarg, t(evalres$rmses), type = "l")
apply(evalres$rmses, 2, max)
