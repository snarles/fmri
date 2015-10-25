####
##  Approximate pchisq(q = r^2, df = d, ncp = theta^2)
####

##source("approximation//polynomial.R")

load("approximation/laplaceData.Rdata")

true_ans <- function(d, theta, r) pchisq(r^2, d, theta^2)
log_true_ans <- function(d, theta, r) pchisq(r^2, d, theta^2, log.p=TRUE)
len_x <- function(theta, r, x) {
  1/theta * sqrt(theta^2 * (theta - r + x)^2
                 - 4 * (theta + x/2) * (theta + x/2  - r) *
                   (r - x/2) * (x/2))
}
ip_x <- function(theta, r, x) len_x(theta, r, x)/(theta - r + x)
q_x <- function(theta, r, x) theta^2 * (theta - r + x)^2 -
  4 * (theta + x/2) * (theta + x/2 - r) * (r - x/2) * x/2
dchi <- function(x, df) 2^(1 - (df/2)) * x^(df - 1) * exp(-x^2/2)/gamma(df/2)
density_x <- function(d, theta, r, x) {
  dchi(theta - r + x, d) * pbeta((1 - ip_x(theta, r, x))/2, (d-1)/2, (d-1)/2)
}
log_density_x <- function(d, theta, r, x) {
  (1 - (d/2)) * log(2) - lgamma(d/2) +
    (d-1) * log(theta - r + x) - (theta - r + x)^2/2 + 
    pbeta((1 - ip_x(theta, r, x))/2, (d-1)/2, (d-1)/2, log.p = TRUE)
}
logsumexp <- function(v) log(sum(exp(v - max(v)))) + max(v)
logsumexp_sign <- function(ss, ls) {
  m <- max(ls)
  s0 <- sum(exp(ls - m) * ss)
  c(s = sign(s0), l = log(abs(s0)))
}
logminusexp <- function(a, b) {
  m <- max(c(a, b))
  s <- sign(a - b)
  ans <- log(abs(exp(a - m) - exp(b - m))) + m  
  c(s = s, l = ans)
}
logmult <- function(v) {
  s <- prod(sign(v))
  c(s = s, l = sum(log(abs(v))))
}

## numerical intergration
prox_ans <- function(d, theta, r, delta = 1e-3) {
  xs <- seq(0, r, by = delta)
  sum(density_x(d, theta, r, xs)) * delta
}
log_prox_ans <- function(d, theta, r, delta = 1e-3) {
  xs <- seq(0, r, by = delta)
  logsumexp(log_density_x(d, theta, r, xs)) + log(delta)
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
  suppressWarnings(
    ans <- 1/beta(a, b) * ((x^(a-2) * (1-x)^(b-1))*(a-1) - 
                    (x^(a-1) * (1-x)^(b-2))*(b-1))
  )
  if (sum(is.na(ans)) > 0) {
    res <- logminusexp((a-2) * log(x) + (b-1)*log(1-x) + log(a-1),
                       (a-1) * log(x) + (b-2)*log(1-x) + log(b-1))
    sn <- res[1]
    mag <- -lbeta(a, b) + res[2]
    ans <- sn * exp(mag)
  }
  ans
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

deriv2_logFp <- function(d, theta, r, x, heavy = TRUE) {
  ff <- dbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
  FF <- pbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
  fp <- deriv_dbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
  di <- deriv_ip_x(theta, r, x)
  di2 <- deriv2_ip_x(theta, r, x)
  if (!heavy)
    return(-ff^2/FF^2 * di^2/4 + fp/FF * di^2/4 - ff/FF * di2/2)
  res1 <- logmult(c(ff, 1/FF, di))
  res1[1] <- 1
  res1[2] <- res1[2]*2 - log(4)
  res2 <- logmult(c(fp, 1/FF, di^2, 1/4))
  res3 <- logmult(c(ff, 1/FF, di2, 1/2))
  res <- logsumexp_sign(c(res1[1], res2[1], res3[1]),
                        c(res1[2], res2[2], res3[2]))
  res[1] * exp(res[2])
}

## laplace approximation
log_prox_ans2 <- function(d, theta, r) {
  objective_f <- function(x) log_density_x(d, theta, r, x)
  res <- optimise(objective_f, interval = c(0, r), maximum = TRUE)
  x_star <- res$maximum
  log_g_star <- res$objective
  omega <- -deriv2_logFp(d, theta, r, x_star) -
    deriv2_log_dchi(theta + r - x_star, d)
  log_g_star + 1/2 * (log(2 * pi) - log(omega))
}

laplace_x_star <- function(d, theta, r) {
  ratio <- (r/theta)^2
  ind <- order(abs(ratio - laplace_ratios))[1]
  laplace_coeffs[ind] * sqrt(d) + laplace_ints[ind]
}

## laplace approximation : plug in the mysterious x* = sqrt(d)
log_prox_ans3 <- function(d, theta, r) {
  x_star <- laplace_x_star(d, theta, r)
  log_g_star <- log_density_x(d, theta, r, x_star)
  omega <- -deriv2_logFp(d, theta, r, x_star) -
    deriv2_log_dchi(theta + r - x_star, d)
  if (is.na(omega)) omega <- 1 ## just ignore it
    #print(paste("d=", d, ";theta=", theta, ";r=", r))
  if (omega < 0) omega <- 1  ## just ignore it
    #print(paste("d=", d, ";theta=", theta, ";r=", r))
  log_g_star + 1/2 * (log(2 * pi) - log(omega))
}

pchisq_laplace <- function(q, df, ncp) {
  theta <- sqrt(ncp)
  r <- sqrt(q)
  ans <- min(c(1, exp(log_prox_ans3(df, theta, r))))
  ans
}

## just use volume approximation!! duh!!
log_prox_ans4 <- function(d, theta, r) {
  ## log-volume of sphere
  lv <- (d/2) * log(2 * pi) - lgamma(d/2) + log(r) - log(d)
  ## log-density of gaussian
  ld <- -d/2 * log(2 * pi) - (theta-r)/2
  lv + ld
}