## Extrapolation method based on Kernel Density estimation

logsumexp <- function(v) max(v) + log(sum(exp(v - max(v))))

## KDE estimate of f(x), then get CDF F(x)
log_gaussian_kernel_cdf <- function(xs, x, bw = "ucv", return.bw = FALSE) {
  dens <- density(xs, bw=bw)
  ps <- pnorm(x, mean = xs, sd = dens$bw, log.p=TRUE)
  if (return.bw) return(dens$bw)
  logsumexp(ps) - log(length(ps))
}

log_raccs <- function(pmat, ...) {
  sapply(1:nrow(pmat), function(ind) log_gaussian_kernel_cdf(pmat[ind, -ind], pmat[ind, ind], ...))
}

kernel_extrap3 <- function(pmat, Ks, ...) {
  log_racs <- log_raccs(pmat, ...)
  sapply(Ks, function(k) mean(exp((k-1)*log_racs)))
}


# ind <- 3
# xs <- pmat[ind, -ind]
# x <- pmat[ind, ind]
# bw = "bcv"


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
