####
## METHOD FOR MUTUAL INFORMATION ESTIMATION
###


####
##  OUR METHOD: I_LOWSNR
####


meanexp <- function(v) {
  vm <- max(v)
  mean(exp(v - vm)) * exp(vm)
}

## function f_K
fK <- function(mus, K, mc.reps = 1e4, naive = FALSE) {
  if (naive) {
    xs <- randn(mc.reps, K)
    xs[, 1] <- xs + mus[1]
    res <- apply(xs, 1, function(v) order(v)[1])
    sum(res != 1)/mc.reps
  }
  samp <- qnorm(((1:mc.reps) - 0.5)/mc.reps)
  sampmat <- repmat(t(samp), length(mus), 1) - mus# one row per mu  
  temp <- log(1 - pnorm(sampmat))
  1 - apply((K-1) * temp, 1, meanexp)
}

## inverse fK
inv_fK <- function(p, K, upper = 10, res = 1e3) {
  if (p ==0) return(Inf)
  xs <- seq(0, upper, length.out = res + 1)
  ps <- fK(xs, K)
  xs[order(abs(ps- p))[1]]
}

Ihat_LS <- function(p, K, upper = 10, res = 1e3) {
  1/2 * inv_fK(p, K)^2
}

## Ihat = 2 inv_fK^2

#inv_fK(0.2, 10)
#fK(2.46, 10)

####
##  FANO METHOD
####

Ihat_fano <- function(p, K) {
  if (p == 0) return(Inf)
  log(K) - p*log(K-1) + p*log(p) + (1-p)*log(1-p)
}

####
##  CM METHOD
####

cm_to_I <- function(cm) {
  m1 <- rowSums(cm)
  m2 <- colSums(cm)
  pm <- m1 %*% t(m2)
  temp <- 0*cm
  temp[cm != 0] <- cm[cm != 0] * log(cm[cm != 0]/pm[cm != 0])
  sum(temp)
}
