library(parallel)
library(abind)

mcc <- 40

source("approximation/gaussian_identity_finsam.R")
source("extrapolation/ku_source.R")

p <- 10


sigma2s <- seq(0.3, 0.7, 0.05)^2
# ac_by_sig2 <- sapply(sigma2s, function(s) {
#   1 - mc_ident_fs(p, s, s, K, 10)
# })
# plot(ac_by_sig2)
# sigma2 <- 0.1
# sigma2_tr <- sigma2 ## equivalent to 1 nn
K <- 2000
ksub <- 500
#mc.reps <- 2
# 
# mus <- randn(K, p)
# ys <- mus + sqrt(sigma2) * randn(K, p)
# mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
# pmat <- -pdist2(mu_hats, ys)
# accs <- 1 - resample_misclassification(pmat, 1:K, 1:K)
# plot(accs, type = "l")

nsplines <- 100
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- rev(1 - knts[-c(1, nsplines + 2)]^2)
MM_2_sq <- spline1_moments(knts, 1:K)

nsplines <- 100
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- knts[-c(1, nsplines + 2)]
MM_2 <- spline1_moments(knts, 1:K)

nsplines <- 1000
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- knts[-c(1, nsplines + 2)]
MM_3 <- spline1_moments(knts, 1:K)

nsplines <- 1000
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- rev(1 - knts[-c(1, nsplines + 2)]^2)
MM_3_sq <- spline1_moments(knts, 1:K)

nsplines <- 10000
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- knts[-c(1, nsplines + 2)]
MM_4 <- spline1_moments(knts, 1:K)

nsplines <- 10000
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- rev(1 - knts[-c(1, nsplines + 2)]^2)
MM_4_sq <- spline1_moments(knts, 1:K)

basis_vecs <- list(mm2 = MM_2, mm2sq =MM_2_sq, mm3 = MM_3, mm3sq = MM_3_sq, mm4 = MM_4, mm4sq = MM_4_sq)
all_final_preds <- array(0, dim = c(mc.reps, length(ksubs), length(basis_vecs)))
all_accs <- matrix(0, mc.reps, K)


subfun <- function(repno) {
  set.seed(repno)
  preds <- matrix(0, length(sigma2s), length(basis_vecs))
  mus <- randn(K, p)
  noise1 <- randn(K, p)
  noise2 <- randn(K, p)
  accsS <- matrix(0, length(sigma2s), K)
  for (j in 1:length(sigma2s)) {
    sigma2 <- sigma2s[j]
    sigma2_tr <- sigma2s[j]
    ys <- mus + sqrt(sigma2) * noise1
    mu_hats <- mus + sqrt(sigma2_tr) * noise2
    pmat <- -pdist2(mu_hats, ys)
    accs <- 1 - resample_misclassification(pmat, 1:K, 1:K)
    accsS[j, ] <- accs
    pmat_sub <- pmat[1:ksub, 1:ksub]
    accs_sub <- 1 - resample_misclassification(pmat_sub, 1:ksub, 1:ksub)
    for (ind in 1:length(basis_vecs)) {
      MM <- basis_vecs[[ind]]
      bt <- nnls::nnls(MM[1:ksub, ], 1- accs_sub)
      preds[j, ind] <- 1 - (MM[K, , drop = FALSE] %*% bt$x)
    }
  }
  list(preds = preds, accsS = accsS)
}

t1 <- proc.time()
res <- mclapply(1:80, subfun, mc.cores = 10)
proc.time() - t1

all_accs <- lapply(res, `[[`, "accsS")
all_accs <- do.call(abind, c(all_accs, list(rev.along = 0)))
dim(all_accs)

all_final_preds <- lapply(res, `[[`, "preds")
all_final_preds <- do.call(abind, c(all_final_preds, list(rev.along = 0)))


true_accs <- apply(all_accs, c(1, 2), mean)

matplot(t(true_accs), type = "l", ylim = c(0, 1))

errs <- (all_final_preds - true_accs[, K])^2

final_pred_err_mat <- apply(errs, c(1, 2), mean)
pred_mat <- apply(all_final_preds, c(1, 2), mean)

matplot(true_accs[, K], pred_mat, type = "l", xlim = c(0, 1), ylim = c(0, 1), lwd = 3)
abline(0, 1, col = "grey", lwd = 3)

#final_pred_err_mat <- apply(errs[, , sample(80, replace = TRUE)], c(1, 2), mean)
matplot(log(final_pred_err_mat), type = "l", lwd = 3)
legend(10, 0, legend = names(basis_vecs), col = 1:6, lty = 1:6)

save(all_accs, errs, true_accs, all_final_preds, final_pred_err_mat, file = "approximation/mcgs2fovs.rds")



