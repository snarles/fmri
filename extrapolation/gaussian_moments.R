require(pracma)

logsumexp <- function(v) log(sum(exp(v - max(v)))) + max(v)
logmeanexp <- function(v) logsumexp(v - log(length(v)))
# log(mean(exp(1:5)))
# logmeanexp(1:5)

# mus <- seq(-2,2,.5)
# sigma2 <- 1
# 
# mu <- 1
# 
# ks <- seq(100, 5000, 100)

gaussian_basis <- function(mu, sigma2, ks, mc.reps = 1e4) {
  samp_qs <- (1:mc.reps - 0.5)/mc.reps
  samp <- qnorm(samp_qs, mu, sqrt(sigma2))
  #c(sd(samp), sqrt(sigma2), sd(samp) - sqrt(sigma2))
  logps <- pnorm(samp, log.p = TRUE)
  Logps <- repmat(logps, length(ks), 1)
  Ks <- repmat(t(t(ks - 1)), 1, mc.reps)
  ans <- apply(Logps * Ks, 1, logmeanexp)
  exp(ans)
}

# mu <- 5
# sigma2 <- 0.5
# ks <- seq(100, 5000, 100)
# ans1 <- gaussian_basis(mu, sigma2, ks, 1e6)
# ans2 <- gaussian_basis(mu, sigma2, ks, 1e5)
# max(abs(ans1 - ans2))
