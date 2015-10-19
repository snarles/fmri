source("approximation/large_deviations_source.R")

sigma2 <- 1
ds <- c(1:10, 100 * 1:30, 1e4 * 1:10)
ress <- matrix(0, length(ds), 3)
for (i in 1:length(ds)) {
  d <- ds[i]
  theta <- sqrt((1 + sigma2) * d); r <- sqrt(sigma2 * d)
  ress[i, ] <-   c(log_true_ans(d, theta, r),
  log_prox_ans(d, theta, r),
  log_prox_ans2(d, theta, r))  
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
      ress[i, j] <- x_star/r
    }
  }
  ress
}

ds <- 100 * 1:100
sigma2s <- seq(0.1, 10, by = 0.1)^2
ress <- ratiodata(sigma2s, ds)
matplot(ress, type = "l")
plot(sigma2s, ress[dim(ress)[1], ], type = "l", ylim = c(0, 1))
plot(sigma2s, 1/ress[dim(ress)[1], ], type = "l")
plot(sigma2s, 1/ress[dim(ress)[1], ]^2, type = "l")
