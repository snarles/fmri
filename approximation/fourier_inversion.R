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

ts <- seq(-10, 10, length.out = 1e4)

dx <- 0.1
xs <- seq(0, 40, by = dx)
betas <- c(2, 1, 1)
samp <- ss_samp(betas)
mean(samp < 40)

fx <- finv(betas, xs)
sum(fx) * dx

mean(pnorm(sqrt(samp)/2))
sum(fx * dx * pnorm(sqrt(xs)/2))
plot(xs, fx, type = "l")

plot(ts, Re(chif(betas, ts)), type = 'l')
plot(ts, Im(chif(betas, ts)), type = 'l')
plot(ts, abs(chif(betas, ts)), type = 'l')
