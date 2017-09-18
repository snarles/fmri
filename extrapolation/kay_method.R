## Extrapolation method based on Kernel Density estimation

## KDE estimate of f(x), then get CDF F(x)
gaussian_kernel_cdf <- function(xs, x, bw = "bcv", return.bw = FALSE) {
  dens <- density(xs, bw=bw)
  ps <- pnorm(x, mean = xs, sd = dens$bw)
  if (return.bw) return(dens$bw)
  mean(ps)
}

raccs <- function(pmat, ...) {
  sapply(1:nrow(pmat), function(ind) gaussian_kernel_cdf(pmat[ind, -ind], pmat[ind, ind], ...))
}

kernel_extrap <- function(pmat, Ks, ...) {
  racs <- raccs(pmat, ...)
  sapply(Ks, function(k) mean(racs^(k-1)))
}


## test it out
# source("approximation/gaussian_identity_finsam2.R")
# p <- 10
# K <- 1e4
# ksub <- 1e3
# kz <- (1:1e2)*1e2
# sigma2 <- 0.3
# sigma2_tr <- sigma2
# mus <- randn(K, p)
# ys <- mus + sqrt(sigma2) * randn(K, p)
# mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
# pmat <- -pdist2(ys, mu_hats)
# rSqs <- rowSums((ys - mu_hats)^2)
# counts <- countDistEx(mu_hats, ys, rSqs)
# accs <- count_acc(counts, kz)    
# accs_kern <- kernel_extrap(pmat, kz)
# plot(kz, accs, type = "l", ylim = c(0, 1))
# lines(kz, accs_kern, col = "red")
