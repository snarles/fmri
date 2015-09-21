## approximating chi-squared sums sum beta_i X_i^2

library(pracma)

ss_samp <- function(betas, n) {
  p <- length(betas)
  colSums(betas * randn(p, n)^2)
}

kfunc0 <- function(betas, a) {
  p <- length(betas)
  amat <- repmat(t(t(a)), 1, p)
  bmat <- repmat(betas, length(a), 1)
  -1/(2 * length(betas)) * rowSums(log(1 - 2*amat*bmat))
}
  
kfunc1 <- function(betas, a) {
  p <- length(betas)
  amat <- repmat(t(t(a)), 1, p)
  bmat <- repmat(betas, length(a), 1)
  1/(2 * p) * rowSums(1/(amat - (1/(2 * bmat))))
}

kfunc2 <- function(betas, a) {
  p <- length(betas)
  amat <- repmat(t(t(a)), 1, p)
  bmat <- repmat(betas, length(a), 1)
  -1/(2 * length(betas)) * rowSums((1/(amat - (1/(2 * bmat))))^2)  
}

build_a_table <- function(betas, lb = 1e-3, ub = 4*sum(betas), gap = 1e-2) {
  mb <- min(1/(2 * betas))
  la <- 1 + mb; ua <- 2 + mb
  while(kfunc1(betas, ua) > lb) ua <- 2 * ua
  while(kfunc1(betas, la) < ub) la <- (la - mb)/2 + mb
  as <- seq(la, ub, length.out = 100)
  tab <- data.frame(as = as, ts = kfunc1(betas, as))
  flag <- TRUE
  while (flag) {
    temp <- tab$ts[-length(tab$ts)] - tab$ts[-1]
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
  tab
}

betas <- c(2,1,1)
tab <- saddleprox(betas)

View(tab)
kfunc2(betas, 3)
