## approximating chi-squared sums sum beta_i X_i^2

library(pracma)

ss_samp <- function(betas, n = 1e5) {
  p <- length(betas)
  colSums(betas * randn(p, n)^2)
}

chif <- function(betas, ts) {
  p <- length(betas)
  amat <- repmat(t(t(ts)), 1, p)
  bmat <- repmat(betas, length(ts), 1)
  cfm <- -1/2 * log(1 - 2 * bmat * amat * (1i))
  exp(rowSums(cfm))
}

finv <- function(betas, xs, ts = seq(-10, 10, length.out = 1e4)) {
  dt <- ts[2] - ts[1]
  cf <- chif(betas, ts)
  tmat <- repmat(ts, length(xs), 1)
  xmat <- repmat(t(t(xs)), 1, length(ts))
  fmat <- exp((-1i) * tmat * xmat)
  as.numeric(Re(dt/(2 * pi) * fmat %*% cf))
}

dx <- 0.1
xs <- seq(0, 20, by = dx)
betas <- c(2, 1, 1)
samp <- ss_samp(betas)
fx <- finv(betas, xs)
sum(fx) * dx
mean(samp < 20)

