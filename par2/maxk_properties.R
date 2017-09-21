## mean and variance of max_k gaussians

require(pracma)

gaussian_maxk_moments <- function(k, naive = FALSE, mc.reps = 1e4) {
  if (naive) {
    zs <- apply(randn(mc.reps, k), 1, max)
    return(
      c(mean(zs), var(zs))
      )
  }
  grid <- seq(0, log(mc.reps), length.out = mc.reps)
  dlt <- grid[2] - grid[1]
  ps_p <- pnorm(grid)
  ps_n <- pnorm(-grid)
  (m1 <- sum((1 - (ps_p^k)) * dlt)  - sum((ps_n^k) * dlt))
  (m2 <- 2 * (sum(grid * (1 - (ps_p^k)) * dlt)  + sum(grid * (ps_n^k) * dlt)))
  (vv <- m2 - m1^2)
  c(m1, vv)
}

approx_gaussian_maxk_moments <- function(k) {
  if (k == 1) return(c(0, 1))
  mu.coefs <- c(-0.676, 1.05)
  var.coefs <- c(0.048984, 1.15798, 0.001205)
  var.coefs_sm <- c(0.990, 0.648, 0.078)
  mu.hat <- mu.coefs[1] + mu.coefs[2] * sqrt(2 * log(k))
  if (k < 100) {
    pres.hat <- var.coefs_sm[1] + var.coefs_sm[2] * log(k) + var.coefs_sm[3] * log(k)^2
  } else {
    pres.hat <- var.coefs[1] + var.coefs[2] * log(k) + var.coefs[3] * log(k)^2
  }
  c(mu.hat, 1/pres.hat)
}

gaussian_maxk_moments(10, TRUE)
gaussian_maxk_moments(10)

gaussian_maxk_moments(20, TRUE)
gaussian_maxk_moments(20)


ks <- (1:1000)^2
t1 <-proc.time()
moms <- sapply(ks, gaussian_maxk_moments, mc.reps = 1e6)
proc.time() - t1




dim(moms)
save(ks, moms, file = "par2/maxk_properties_precomp.RData")
load("par2/maxk_properties_precomp.RData")


plot(sqrt(2 * log(ks)), moms[1, ], type = "l")
plot(log(ks), 1/moms[2, ], type = "l")

filt <- 10:length(ks)
summary(lm(moms[1, filt] ~ sqrt(2 * log(ks[filt]))))


mu.coefs <- c(-0.676, 1.05)

plot(log(ks[filt]), 1/moms[2, filt], type = "l")

summary(lm(1/moms[2, filt] ~ log(ks[filt]) + I(log(ks[filt])^2)))

coef(lm(1/moms[2, filt] ~ log(ks[filt]) + I(log(ks[filt])^2)))

var.coefs <- c(-0.0677958, 1.1824989)
var.coefs <- c(0.048984, 1.15798, 0.001205)

plot(log(ks[filt]), (1/moms[2, filt])-(var.coefs[1] + var.coefs[2] * log(ks[filt])), type = "o")
plot(log(ks[filt]), (1/moms[2, filt])-(var.coefs[1] + var.coefs[2] * log(ks[filt]) + var.coefs[3] * log(ks[filt])^2), type = "o")



## find transition point

plot(diff(moms[1, ])/diff(sqrt(2 * log(ks))))
plot(diff(1/moms[2, ])/diff(log(ks)))




ks <- (1:200)
moms <- sapply(ks, gaussian_maxk_moments, mc.reps = 1e5)
apmoms <- sapply(ks, approx_gaussian_maxk_moments)
plot(moms[1, ], apmoms[1, ], type = "o")
abline(0, 1, col = "red")

plot(moms[2, ], apmoms[2, ], type = "o")
abline(0, 1, col = "red")


plot(sqrt(2 * log(ks)), moms[1, ], type = "o")
plot(log(ks), 1/moms[2, ], type = "o")
summary(lm(1/moms[2, ] ~ log(ks) + I(log(ks)^2)))

var.coefs_sm <- c(0.834, 0.912)
var.coefs_sm <- c(0.990, 0.648, 0.078)


plot(sqrt(2 * log(ks)), moms[1, ] - mu.coefs[1] - mu.coefs[2] * sqrt(2 * log(ks)), type = "o")
plot(log(ks), (1/moms[2, ])-(var.coefs_sm[1] + var.coefs_sm[2] * log(ks)), type = "o")
plot(log(ks), (1/moms[2, ])-(var.coefs_sm[1] + var.coefs_sm[2] * log(ks) + var.coefs_sm[3] * log(ks)^2), type = "o")


plot(log(ks[-1]), (1/moms[2, -1])/(var.coefs[1] + var.coefs[2] * log(ks[-1])), type = "o")


plot(diff(moms[1, ])/diff(sqrt(2 * log(ks))), type = "o")
plot(diff(1/moms[2, ])/diff(log(ks)), type = "o")
