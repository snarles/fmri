library(lineId)
library(pracma)

binmom <- function(succ, tot, k) {
  choose(succ, k)/choose(tot, k)
}

get_avrisk_mat <- function(ks, d) {
  ans <- matrix(0, length(ks), d + 1)
  ans <- col(ans) - 1
  ans <- ans + (ks - 1)
  ans <- 1/ans
  ans <- ans * (ks - 1)
  ans
}

get_kmat <- function(k, d) {
  kmat <- matrix(0, k - 1, d + 1)
  al <- row(kmat)
  bt <- k - al
  h <- col(kmat) - 1
  kmat <- lgamma(al + h) + lgamma(al + bt)  - lgamma(al) - lgamma(al + bt + h)
  kmat <- exp(kmat)
  kmat  
}

get_rank_prop <- function(pmat, true_ys) {
  k <- nrow(pmat)
  p2 <- apply(pmat, 2, rank)
  true_ranks <- p2[cbind(true_ys, 1:ncol(pmat))]
  tab <- table(true_ranks)
  ans <- numeric(k)
  ans[as.numeric(names(tab))] <- tab
  ans <- ans/ncol(pmat)
  cumsum(ans)[1:(k-1)]
}

get_sub_errs <- function(pmat, true_ys, ks) {
  k <- nrow(pmat)
  p2 <- apply(pmat, 2, rank)
  true_ranks <- p2[cbind(true_ys, 1:ncol(pmat))]
  1 - sapply(ks, function(v) mean(binmom(true_ranks - 1, k - 1, v - 1)))
}

get_vande <- function(xs = seq(0, 1, 0.01), d) {
  ans <- xs %*% t(rep(1, d+1)) 
  ans <- ans ^ (col(ans) - 1)
  ans
}

avrisk <- function(k, bt) {
  d <- length(bt) - 1
  sapply(k, function(k) (k - 1) * sum(bt/(0:d + k - 1)))
}

## spline which is [t-t0]_+
spline1_maker <- function(t0) {
  function(x) pmax(0, x - t0)
}

get_f_moments <- function(ff, ks, res = 1e3) {
  ts <- seq(0, 1, 1/res)
  ys <- ff(ts)
  Tmat <- sapply(ks, function(k) (k - 2) * ts^(k - 2))
  Tprod <- Tmat * ys
  fmoms <- colSums(Tprod)/res
  fmoms
}

make_moment_mat <- function(flist, ks, res = 1e3) {
  sapply(flist, function(fl) get_f_moments(fl, ks, res))
}

# ffs <- lapply(seq(0, 1, 0.1), spline1_maker)
# plot(ffs[[5]](seq(0, 1, 0.01)))
# 
# ks <- 4:10
# make_moment_mat(ffs, ks, res = 1e5)
