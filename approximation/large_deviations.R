####
##  Approximate pchisq(q = r^2, df = d, ncp = theta^2)
####

true_ans <- function(d, theta, r) pchisq(r^2, d, theta^2)
len_x <- function(theta, r, x) {
  1/theta * sqrt(theta^2 * (theta - r + x)^2
                  - 4 * (theta + x/2) * (theta + x/2  - r) *
                    (r - x/2) * (x/2))
}
ip_x <- function(theta, r, x) len_x(theta, r, x)/(theta - r + x)
dchi <- function(x, df) 2^(1 - (df/2)) * x^(df - 1) * exp(-x^2/2)/gamma(df/2)
density_x <- function(d, theta, r, x) {
  dchi(theta - r + x, d) * pbeta((1 - ip_x(theta, r, x))/2, (d-1)/2, (d-1)/2)
}
log_density_x <- function(f, theta, r, x) {
  (1 - (df/2)) * log(2) - lgamma(df/2) +
    (d-1) * log(theta - r + x) - (theta - r + x)^2/2 + 
    pbeta((1 - ip_x(theta, r, x))/2, (d-1)/2, (d-1)/2, log.p = TRUE)
}
## numerical intergration
prox_ans <- function(d, theta, r, delta = 1e-3) {
  xs <- seq(0, r, by = delta)
  sum(density_x(d, theta, r, xs)) * delta
}
## laplace approximation
prox_ans2 <- function(d, theta, r) {
  objective_f <- function(x) log(density_x(d, theta, r, x))
  res <- optimise(objective_f, interval = c(0, r), maximum = TRUE)
  x_star <- res$maximum
  log_g_star <- res$objective
}


## check len function
d <- 10; theta <- 5; r<- 1
x <- 0
(len <- len_x(theta, r, x))
ip_x(theta, r, x)
y <- sqrt((theta - r + x)^2 - len^2)
c(sqrt(r^2 - y^2) + len, theta)

## approximation error
d <- 20; theta <- 40; r<- 10
true_ans(d, theta, r)
prox_ans(d, theta, r)

## plot density x
delta <- 1e-2
xs <- seq(0, r, by = delta)
plot(xs, ip_x(theta, r, xs), type = 'l')
plot(xs, density_x(d, theta, r, xs), type = "l")

## plot log-density
d <- 80; theta <- 200; r<- 1
delta <- 1e-2
xs <- seq(0, r, by = delta)
ip_x(theta, r, xs)
ds <- log_density_x(d, theta, r, xs)
md <- max(ds)
plot(xs[ds > md - 10], ds[ds > md - 10], type = "l")


