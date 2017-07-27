library(parallel)
library(abind)

mcc <- 40

source("approximation/gaussian_identity_finsam.R")
source("extrapolation/ku_source.R")

p <- 10



p <- 10
sigma2 <- 0.3
sigma2_tr <- sigma2 ## equivalent to 1 nn
K <- 2000
ksub <- 500

load("approximation/mcgs2f.rda", verbose = TRUE) ## get the true_acces
matplot(t(all_accs), type = "l")
ses <- apply(all_accs, 2, sd)/sqrt(nrow(all_accs))
plot(1:K, true_accs, type = "l", ylim = c(0, 1))
lines(1:K, true_accs + ses, col = "red")
lines(1:K, true_accs - ses, col = "red")


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

t1 <- proc.time()
repno <- 1
set.seed(repno)
mus <- randn(K, p)
noise1 <- randn(K, p)
noise2 <- randn(K, p)
ys <- mus + sqrt(sigma2) * noise1
mu_hats <- mus + sqrt(sigma2_tr) * noise2
pmat <- -pdist2(mu_hats, ys)
accs <- 1 - resample_misclassification(pmat, 1:K, 1:K)
pmat_sub <- pmat[1:ksub, 1:ksub]
accs_sub <- 1 - resample_misclassification(pmat_sub, 1:ksub, 1:ksub)
preds <- matrix(0, length(basis_vecs), K)
for (ind in 1:length(basis_vecs)) {
  MM <- basis_vecs[[ind]]
  bt <- nnls::nnls(MM[1:ksub, ], 1- accs_sub)
  preds[ind, ] <- 1 - (MM %*% bt$x)
}
proc.time() - t1


save(true_accs, preds, file = "approximation/mcgs2_oneplot.rda")

pdf("approximation/mcgs2_oneplot.pdf", height = 4, width = 6)
orig_cols <- hsv(1:3/3,
                 rep(1, 3),
                 0.3 + 0.7 * (1:3/3))
cols <- rep(orig_cols, each = 2)
ltys <- rep(c(2, 3), 3)

plot(1:K, true_accs, type = 'l', lwd = 3, xlab = "k", ylab = "AGA", ylim = c(0, 1))
matplot(ksub:K, t(preds[, ksub:K]), type = "l", col = cols, lty = ltys, lwd = 2, add = TRUE)
abline(v = ksub, lty = 2)
legend(1000, 1, c("true", names(basis_vecs)), col = c("black", cols), lty = c(1, ltys), lwd = 2)
dev.off()
