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


ratiodata <- function(sigma2, ds) {
  rs <- ds * 0; xs <- ds * 0
  for (i in 1:length(ds)) {
    d <- ds[i]
    theta <- sqrt((1 + sigma2) * d); r <- sqrt(sigma2 * d)
    objective_f <- function(x) log_density_x(d, theta, r, x)
    res <- optimise(objective_f, interval = c(0, r), maximum = TRUE)
    x_star <- res$maximum
    rs[i] <- r; xs[i] <- x_star
  }
  data.frame(ds, rs, xs, ratio = xs/rs)  
}

sigma2 <- 1
ds <- c(1:10, 100 * 1:20)

View(ratiodata(1, ds))
View(ratiodata(2, ds))
View(ratiodata(3, ds))
