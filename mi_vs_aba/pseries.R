## computation of I vs ABA curves using power series

logsumexp <- function(x) max(x) + log(sum(exp(x - max(x))))

cgf_k <- function(bts, k, nterms = 200) {
  bts <- t(bts)
  js <- 0:nterms
  Bm <- rep(1, length(js)) %*% bts
  Jm <- js %*% t(rep(1, length(bts)))
  rawlogterms <- Jm * log(Bm) - lgamma(Jm + 1) - log(Jm + (1/(k-1)))
  apply(rawlogterms, 2, logsumexp) - log(k-1)
}

d_cgf_k <- function(bts, k, nterms = 200) {
  a0 <- cgf_k(bts, k, nterms)
  bts <- t(bts)
  js <- 1:nterms
  Bm <- rep(1, length(js)) %*% bts
  Jm <- js %*% t(rep(1, length(bts)))
  rawlogterms <- (Jm-1) * log(Bm) - lgamma(Jm + 1) - log(Jm + (1/(k-1))) + log(Jm)
  exp(apply(rawlogterms, 2, logsumexp) - a0 - log(k-1))
}

## check derivative
# k <- 5
# ff <- function(b) cgf_k(b, k)
# bt <- rexp(1)
# c(bt, numDeriv::grad(ff, bt), d_cgf_k(bt, k))

## comparison
# bts <- exp(-2:5)
# k <- 3
# 
# source("mi_vs_aba/alpha_vs_beta.R")
# pts2 <- sapply(bts, function(x) {
#   qs <- qs_par(x, k, reso)
#   get_q_ints(qs, k)
# })
# 
# cbind(pts2[3, ], d_cgf_k(bts, k))
# cbind(pts2[2, ], bts * d_cgf_k(bts, k) - cgf_k(bts, k))

