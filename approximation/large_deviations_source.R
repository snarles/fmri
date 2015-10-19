####
##  Approximate pchisq(q = r^2, df = d, ncp = theta^2)
####

source("approximation//polynomial.R")

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
log_density_x <- function(d, theta, r, x) {
  (1 - (d/2)) * log(2) - lgamma(d/2) +
    (d-1) * log(theta - r + x) - (theta - r + x)^2/2 + 
    pbeta((1 - ip_x(theta, r, x))/2, (d-1)/2, (d-1)/2, log.p = TRUE)
}
## numerical intergration
prox_ans <- function(d, theta, r, delta = 1e-3) {
  xs <- seq(0, r, by = delta)
  sum(density_x(d, theta, r, xs)) * delta
}

deriv_log_dchi <- function(x, d) (d - 1)/x - x
deriv2_log_dchi <- function(x, d) -(d - 1)/x^2 - 1
deriv_q_x <- function(theta, r, x) {
  1*x^3 + -3*x^2*r + 2*r^2*x + 3*x^2*theta + 
    -6*theta*r*x + 2*r^2*theta + 4*theta^2*x +
    -4*theta^2*r + 2*theta^3
}
deriv2_q_x <- function(theta, r, x) {
  3*x^2 + -6*r*x + 2*r^2 + 6*theta*x + -6*theta*r + 4*theta^2
}
deriv_sqrt_q <- function(theta, r, x)
  deriv_q_x(theta, r, x)/(2 * sqrt(q_x(theta, r, x)))
deriv2_sqrt_q <- function(theta, r, x) {
  qq <- q_x(theta, r, x)
  -1/4 * deriv_q_x(theta, r, x)^2 * qq^(-3/2) + 
    1/(2 * sqrt(qq)) * deriv2_q_x(theta, r, x)  
}

deriv_dbeta <- function(x, a, b) {
  1/beta(a, b) * ((x^(a-2) * (1-x)^(b-1))*(a-1) - 
                    (x^(a-1) * (1-x)^(b-2))*(b-1))
}
deriv_ip_x <- function(theta, r, x) {
  ip <- ip_x(theta, r, x)
  qq <- q_x(theta, r, x)
  1/theta * (-sqrt(qq)/(theta-r+x)^2 + 
               deriv_sqrt_q(theta, r, x)/(theta-r+x))
}
deriv2_ip_x <- function(theta, r, x) {
  a <- theta - r + x
  1/theta *
    (2 * sqrt(q_x(theta, r, x))/a^3 - 
       2 * deriv_sqrt_q(theta, r, x)/a^2 + 
       deriv2_sqrt_q(theta, r, x)/a)
}
logFp <- function(d, theta, r, x) {
  pbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2, log.p=TRUE)
}

deriv_logFp <- function(d, theta, r, x) {
  ff <- dbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
  FF <- pbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
  di <- deriv_ip_x(theta, r, x)
  -1/2 * ff/FF *di
}

deriv2_logFp <- function(d, theta, r, x) {
  ff <- dbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
  FF <- pbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
  fp <- deriv_dbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
  di <- deriv_ip_x(theta, r, x)
  di2 <- deriv2_ip_x(theta, r, x)
  -ff^2/FF^2 * di^2/4 + fp/FF * di^2/4 - ff/FF * di2/2
}

## laplace approximation
prox_ans2 <- function(d, theta, r) {
  objective_f <- function(x) log(density_x(d, theta, r, x))
  res <- optimise(objective_f, interval = c(0, r), maximum = TRUE)
  x_star <- res$maximum
  log_g_star <- res$objective
  omega <- -deriv2_logFp(d, theta, r, x_star) -
    deriv2_log_dchi(theta + r - x_star, d)
  log_g_star + 1/2 * (log(2 * pi) - log(omega))