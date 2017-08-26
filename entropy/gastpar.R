####
##  Mutual information estimation
####

kldiv <- function(p, q) {
  temp <- p * log(p/q)
  sum(temp[!is.na(temp)])
}

h_mle <- function(freqs) {
  ps <- freqs/sum(freqs)
  temp <- -ps * log(ps)
  sum(temp[!is.na(temp)])
}

mi_naive <- function(ys, h_est = h_mle) {
  k <- dim(ys)[1]
  yl <- sort(unique(as.numeric(ys)))
  tab <- matrix(0, k, length(yl))
  for (i in 1:length(yl)) {
    yy <- yl[i]
    cts <- apply(ys, 1, function(v) sum(v == yy))
    tab[, i] <- cts
  }
  hy <- h_est(colSums(tab))
  ces <- apply(tab, 1, h_est)
  hy - mean(ces)
}

anthropic_correction <- function(ys, alpha) {
  k <- dim(ys)[1]
  yl <- sort(unique(as.numeric(ys)))
  tab <- matrix(0, k, length(yl))
  for (i in 1:length(yl)) {
    yy <- yl[i]
    cts <- apply(ys, 1, function(v) sum(v == yy))
    tab[, i] <- cts
  }
  tT <- colSums(tab)
  pmat <- t(apply(tab, 1, function(v) v/sum(v)))
  ihat <- 0
  for (i in 1:k) {
    pkay <- colSums(((1-alpha)/k + alpha/(k-1)) * pmat) - (alpha)/(k-1) * pmat[i, ]
    ihat <- ihat + kldiv(pmat[i, ], pkay)/k
  }
  ihat  
}

####
##  Data model
####

# k.total <- 200
# ss <- 100
# pp <- matrix(rbeta(k.total * ss, 0.01, 1), k.total, ss)
# pp <- pp/sum(pp)
# 
# px <- rowSums(pp)
# py <- colSums(pp)
# 
# pind <- t(t(px)) %*% t(py)
# 
# (mi_true <- sum(pp * log(pp/pind)))
# 
# k.obs <- 5
# r.each <- 10
# 
# xs <- sample(1:k.total, k.obs, replace = TRUE, prob = px)
# ys <- matrix(0, k.obs, r.each)
# for (i in 1:k.obs) {
#   v <- pp[xs[i], ]; v <- v/sum(v)
#   ys[i, ] <- sample(1:ss, r.each, TRUE, v)
# }
# 
# (mi_0 <- mi_naive(ys))
# (mi_0 <- anthropic_correction(ys, 0))
# (mi_5 <- anthropic_correction(ys, 0.5))
# (mi_9 <- anthropic_correction(ys, 0.9))
# 
# c(mi_true, mi_0, mi_5, mi_9)
