library(parallel)
mcc <- 40

source("approximation/gaussian_identity_finsam.R")
source("extrapolation/ku_source.R")

p <- 10
sigma2 <- 0.3
sigma2_tr <- sigma2 ## equivalent to 1 nn
K <- 2000
ksubs <- 50 * 1:40
mc.reps <- 50

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

basis_vecs <- list(lin2 = MM_2, 
                   lin2q =MM_2_sq, 
                   lin3 = MM_3, 
                   lin3q = MM_3_sq, 
                   lin4 = MM_4, 
                   lin4q = MM_4_sq)
all_final_preds <- array(0, dim = c(mc.reps, length(ksubs), length(basis_vecs)))
all_accs <- matrix(0, mc.reps, K)

t1 <- proc.time()

for (repno in 1:mc.reps) {
  set.seed(repno)
  mus <- randn(K, p)
  ys <- mus + sqrt(sigma2) * randn(K, p)
  mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
  pmat <- -pdist2(mu_hats, ys)
  accs <- 1 - resample_misclassification(pmat, 1:K, 1:K)
  all_accs[repno, ] <- accs
  
  subfun <- function(j) {
    ksub <- ksubs[j]
    pmat_sub <- pmat[1:ksub, 1:ksub]
    accs_sub <- 1 - resample_misclassification(pmat_sub, 1:ksub, 1:ksub)
    preds <- numeric(length(basis_vecs))
    for (ind in 1:length(basis_vecs)) {
      MM <- basis_vecs[[ind]]
      bt <- nnls::nnls(MM[1:ksub, ], 1- accs_sub)
      preds[ind] <- 1 - (MM[K, , drop = FALSE] %*% bt$x)
    }
    preds
  }
  
  # for (j in 1:length(ksubs)) {
  #   ksub <- ksubs[j]
  #   pmat_sub <- pmat[1:ksub, 1:ksub]
  #   accs_sub <- 1 - resample_misclassification(pmat_sub, 1:ksub, 1:ksub)
  #   for (ind in 1:length(basis_vecs)) {
  #     MM <- basis_vecs[[ind]]
  #     bt <- nnls::nnls(MM[1:ksub, ], 1- accs_sub)
  #     pred <- 1 - (MM[K, , drop = FALSE] %*% bt$x)
  #     all_final_preds[repno, j, ind] <- pred
  #   }
  # }
  
  preds <- mclapply(1:length(ksubs), subfun, mc.cores = mcc)
  preds2 <- do.call(rbind, preds)
  all_final_preds[repno, , ] <- preds2
}

proc.time() - t1

true_accs <- colMeans(all_accs)
errs <- (all_final_preds - true_accs[K])^2

final_pred_err_mat <- apply(errs, c(2, 3), mean)
#final_pred_err_mat <- apply(errs[sample(50, replace = TRUE), , ], c(2, 3), mean)
#load("approximation/mcgs2f.rda")
pdf("approximation/fig_mcgs2f.pdf", width = 5, height = 4)
source("approximation/mcgs2_colscheme.R")
matplot(ksubs, sqrt(final_pred_err_mat), 
        type = "l", ylab = "RMSE", xlab = expression(k[1]), xlim = c(250, 2000),
        ylim = c(0, 0.5), col = cols, lty = ltys, lwd = 3)
legend(1570, 0.53, legend = nms, col = cols, lty = ltys, lwd = 3)
dev.off()

#save(all_accs, errs, true_accs, final_pred_err_mat, file = "approximation/mcgs2f.rda")



