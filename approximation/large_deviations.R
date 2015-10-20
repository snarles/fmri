source("approximation/large_deviations_source.R")

sigma2 <- 1
ds <- c(1:10, 100 * 1:30, 1e4 * 1:10)
ress <- matrix(0, length(ds), 4)
for (i in 1:length(ds)) {
  d <- ds[i]
  theta <- sqrt((1 + sigma2) * d); r <- sqrt(sigma2 * d)
  ress[i, ] <-   c(log_true_ans(d, theta, r),
  log_prox_ans(d, theta, r),
  log_prox_ans2(d, theta, r),
  log_prox_ans3(d, theta, r))  
}

View(cbind(ds, ress))

####
##  Relative location of the max
####


ratiodata <- function(sigma2s, ds) {
  ress <- matrix(0, length(ds), length(sigma2s))
  for (i in 1:length(ds)) {
    for (j in 1:length(sigma2s)) {
      d <- ds[i]; sigma2 <- sigma2s[j]
      theta <- sqrt((1 + sigma2) * d); r <- sqrt(sigma2 * d)
      objective_f <- function(x) log_density_x(d, theta, r, x)
      res <- optimise(objective_f, interval = c(0, r), maximum = TRUE)
      x_star <- res$maximum
      #ress[i, j] <- x_star/r
      ress[i, j] <- x_star
    }
  }
  ress
}

## Code to obtain laplaceData.Rdata
# ds <- 25 * 1:400
# t1 <- proc.time()
# sigma2s <- seq(0.001, 20, by = 0.001)^2
# ress <- ratiodata(sigma2s, ds)
# save(list=c("ress", "sigma2s", "ds"), file="approximation/ratiodata.Rdata")
# proc.time() - t1
# 
# matplot(sqrt(ds), ress, type = "l", col = rainbow(length(sigma2s)))
# plot(sigma2s, ress[ds == 10000, ], type = "l")
# plot(sigma2s/(1 + sigma2s), ress[ds == 10000, ], type = "l")
# 
# lr <- lapply(sigma2s, function(s) {
#   lm(ress[, sigma2s ==s] ~ sqrt(ds))
# })
# laplace_coeffs <- sapply(lr, function(r) r$coefficients[2])
# laplace_ints <- sapply(lr, function(r) r$coefficients[1])
# laplace_ratios <- sigma2s/(1 + sigma2s)
# save(list = c("laplace_coeffs", "laplace_ints", "laplace_ratios"),
#      file = "approximation/laplaceData.Rdata")
# 
# plot(sigma2s, coeffs, type = "l")
# 
# plot(sigma2s, coeffs, type = "l")
# 
# plot(sigma2s, ints, type = "l")
# names(coeffs) <- NULL
# coeffs


## conclusion: x* ~= func(sigma^2) * sqrt(d) + func(sigma^2)


####
##  Size of second derivative can be ignored... true??
##  Looks like size of log(second derivative) is constant in sigma2
####

log_g_vs_omega <- function(d, sigma2) {
  theta <- sqrt((1 + sigma2) * d); r <- sqrt(sigma2 * d)
  x_star <- laplace_x_star(d, theta, r)
  log_g_star <- log_density_x(d, theta, r, x_star)
  omega <- -deriv2_logFp(d, theta, r, x_star) -
    deriv2_log_dchi(theta + r - x_star, d)
  c(log_g_star, log(omega))
}

log_g_vs_omega(100, 1)
log_g_vs_omega(1000, 1)
log_g_vs_omega(100, 100)
log_g_vs_omega(1000, 100)
log_g_vs_omega(10000, 100)
