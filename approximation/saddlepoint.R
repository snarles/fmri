## approximating chi-squared sums sum beta_i X_i^2

library(pracma)

ss_samp <- function(betas, n = 1e5) {
  p <- length(betas)
  colSums(betas * randn(p, n)^2)
}

kfunc0 <- function(betas, a) {
  p <- length(betas)
  amat <- repmat(t(t(a)), 1, p)
  bmat <- repmat(betas, length(a), 1)
  -1/(2 * p) * rowSums(log(1 - 2*amat*bmat))
}
  
kfunc1 <- function(betas, a) {
  p <- length(betas)
  amat <- repmat(t(t(a)), 1, p)
  bmat <- repmat(betas, length(a), 1)
  -1/(2 * p) * rowSums(1/(amat - (1/(2 * bmat))))
}

kfunc2 <- function(betas, a) {
  p <- length(betas)
  amat <- repmat(t(t(a)), 1, p)
  bmat <- repmat(betas, length(a), 1)
  1/(2 * length(betas)) * rowSums((1/(amat - (1/(2 * bmat))))^2)  
}

build_a_table <- function(betas, lb = 1e-3, ub = 4*sum(betas), gap = 1e-2) {
  mb <- min(1/(2 * betas))
  la <- -1; ua <- mb - 1e-1
  while(kfunc1(betas, ua) < ub) ua <- mb - (mb - ua)/2
  while(kfunc1(betas, la) > lb) la <- la * 2
  as <- seq(la, ua, length.out = 100)
  tab <- data.frame(as = as, ts = kfunc1(betas, as))
  flag <- TRUE
  while (flag) {
    temp <- -tab$ts[-length(tab$ts)] + tab$ts[-1]
    inds <- which(temp > gap)
    if (length(inds) == 0) {
      flag <- FALSE
      return(tab)
    }
    as <- (tab$as[inds] + tab$as[inds + 1])/2
    ntab <- data.frame(as = as, ts = kfunc1(betas, as))
    tab <- rbind(tab, ntab)
    tab <- tab[order(tab$as), ]
  }
  tab
}

saddleprox <- function(betas, lb = 1e-3, ub = 4*sum(betas), gap = 1e-2) {
  p <- length(betas)
  tab <- build_a_table(betas, lb, ub, gap)
  as <- tab$as; ts <- tab$ts
  tab$gt <- sqrt(p/(2 * pi * kfunc2(betas, as))) * exp(p * (kfunc0(betas, as) - as * ts))
  ## compute normalization
  min(tab$ts)
  max(tab$ts)
  ts2 <- c(0, tab$ts)
  delts <- ts2[-1] - ts2[-length(ts2)]
  nmlz <- sum(delts * tab$gt)
  tab$gt <- tab$gt/nmlz
  tab
}

## expected value
exval <- function(tab, ff) {
  ts2 <- c(0, tab$ts)
  delts <- ts2[-1] - ts2[-length(ts2)]
  sum(tab$gt * delts * ff(tab$ts))  
}

betas <- c(5,2,1)
tab <- saddleprox(betas)

kfunc2(betas, 3)

ys <- ss_samp(betas)
mean(pnorm(sqrt(ys)/2))
mean(pnorm(.5 * sqrt(ys)/2))
mean((ys < 2))
mean(ys)

exval(tab, function(x) {pnorm(sqrt(x)/2)})
exval(tab, function(x) {pnorm(.5 * sqrt(x)/2)})
exval(tab, function(x) {x < 2})
exval(tab, function(x) x)
