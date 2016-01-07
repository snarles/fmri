library(pracma)

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

cols <- hsv(h = 0:3/4, s = 0.9, v = 0.5)
xs <- seq(-5, 10, by = 0.1)
plot(xs, fK(xs, 2), type = "l", col = cols[1], lwd = 3, xlab = expression(mu), ylab = expression(pi[K]))
abline(v = 0, lwd = 2)
lines(xs, fK(xs, 9), type = "l", col = cols[2], lwd = 3)
lines(xs, fK(xs, 99), type = "l", col = cols[3], lwd = 3)
lines(xs, fK(xs, 999), type = "l", col = cols[4], lwd = 3)
