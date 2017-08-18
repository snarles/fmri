library(parallel)
mcc <- 20

source("approximation/gaussian_identity_finsam.R")
source("extrapolation/ku_source.R")

p <- 10
sigma2 <- 0.3
sigma2_tr <- sigma2 ## equivalent to 1 nn
K <- 2000
ksub <- 1000
mc.reps <- 180

# mus <- randn(K, p)
# ys <- mus + sqrt(sigma2) * randn(K, p)
# mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
# pmat <- -pdist2(mu_hats, ys)
# accs <- 1 - resample_misclassification(pmat, 1:K, 1:K)
# plot(accs, type = "l")

nsplines <- c(100, 400)
nrows <- c(250, 500, 1000)
combmat <- cbind(nsplines = rep(nsplines, each = length(nrows)), 
                 nrows = rep(nrows, length(nsplines)))

basis_vecs <- list()

kseq <- function(nr) {
  interv <- floor(K/nr)
  seq(interv, K, by = interv)
}

for (ii in 1:nrow(combmat)) {
  ns <- combmat[ii, 1]
  nr <- combmat[ii, 2]
  nm <- paste0("s",ns,"r",nr)
  knts <- seq(0, 1, length.out = ns + 2)
  knts <- rev(1 - knts[-c(1, ns + 2)]^2)  
  Xmat <- spline1_moments(knts, kseq(nr))
  Xpred <- spline1_moments(knts, K)
  basis_vecs[[nm]] <- list(Xmat = Xmat, Xpred = Xpred, kseq = kseq(nr))
}

all_final_preds <- array(0, dim = c(mc.reps, length(basis_vecs)))
all_accs <- matrix(0, mc.reps, K)

t1 <- proc.time()

subfun <- function (repno) {
  set.seed(repno)
  mus <- randn(K, p)
  ys <- mus + sqrt(sigma2) * randn(K, p)
  mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
  pmat <- -pdist2(mu_hats, ys)
  accs <- 1 - resample_misclassification(pmat, 1:K, 1:K)
  all_accs[repno, ] <- accs  
  pmat_sub <- pmat[1:ksub, 1:ksub]
  accs_sub <- 1 - resample_misclassification(pmat_sub, 1:ksub, 1:ksub)
  preds <- numeric(length(basis_vecs))
  for (ind in 1:length(basis_vecs)) {
    Xmat <- basis_vecs[[ind]]$Xmat
    ks <- basis_vecs[[ind]]$ks
    Xpred <- basis_vecs[[ind]]$Xpred
    bt <- nnls::nnls(Xmat, 1- accs_sub[ks])
    preds[ind] <- 1 - (Xpred %*% bt$x)
  }
  preds
}  

## !!TBC!!

t1 <- proc.time()
preds <- mclapply(1:mc.reps, subfun, mc.cores = mcc)
preds2 <- do.call(rbind, preds)
all_final_preds[repno, , ] <- preds2
proc.time() - t1
# "12.5 per hour"

true_accs <- colMeans(all_accs)
errs <- (all_final_preds - true_accs[K])^2

final_pred_err_mat <- apply(errs, c(2, 3), mean)
#final_pred_err_mat <- apply(errs, c(2, 3), median)
#plot(all_final_preds[,40,5])
#final_pred_err_mat <- apply(errs[sample(mc.reps, replace = TRUE), , ], c(2, 3), mean)
#load("approximation/mcgs2f.rda")
pdf("approximation/fig_mcgs2f.pdf", width = 5, height = 4)
source("approximation/mcgs2_colscheme.R")
matplot(ksubs, sqrt(final_pred_err_mat), 
        type = "l", ylab = "RMSE", xlab = expression(k[1]), xlim = c(250, 2000),
        ylim = c(0, 0.3), col = cols, lty = ltys, lwd = 3)
# matplot(ksubs, sqrt(final_pred_err_mat), 
#         type = "l", ylab = "RMSE", xlab = expression(k[1]),
#         col = cols, lty = ltys, lwd = 3, ylim = c(0, 0.5))

legend(1570, 0.32, legend = nms, col = cols, lty = ltys, lwd = 3)
dev.off()

#save(all_accs, errs, true_accs, final_pred_err_mat, file = "approximation/mcgs2f.rda")



