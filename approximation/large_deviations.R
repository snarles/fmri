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
  ## TO BE CONTD
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

####
##  Check overall approximation
####



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

####
##  Check derivatives of chi
####

d <- 20 * runif(1); xs <- seq(0, d, by = 0.1)
x <- d * runif(1)

plot(xs, log(dchi(xs, d)), type = "l", lwd = 2)
lines(xs, log(dchi(x, d)) + (xs - x) *
     deriv_log_dchi(x, d), col = "red")
lines(xs, log(dchi(x, d)) + (xs - x) *
        deriv_log_dchi(x, d) + 
        (xs - x)^2/2 * deriv2_log_dchi(x, d), col = "blue")

####
##  Check derivatives of q_x
####

q_x <- function(theta, r, x) theta^2 * (theta - r + x)^2 -
  4 * (theta + x/2) * (theta + x/2 - r) * (r - x/2) * x/2
# ip_x = ip_x_2
#theta <- 10; r <- 1; xs <- seq(0, r, by = 1e-2)
#cbind(ip_x(theta, r, xs), ip_x_2(theta, r, xs))
ip_x_2 <- function(theta, r, x) sqrt(q_x(theta, r, x))/theta/(theta - r + x)
# to be simplified
deriv_q_x_raw <- function(theta, r, x)
  1e10 * (q_x(theta, r, x+1e-10) - q_x(theta, r, x))


## Simplify q_x
vv <- AlgDesign::gen.factorial(7, 3)
colnames(vv) <- c("theta", "r", "x")
V <- form_poly_matrix(vv, 4)
vals <- apply(vv, 1, function(v) q_x(v[1], v[2], v[3]))
res <- lm(vals ~ 0+ V)
summary(res)
check_divis(res$coefficients, 4)
display_div(res$coefficients, 4)
# vvPx4       vvPr1x3       vvPr2x2   vvPtheta1x3 vvPtheta1r1x2 vvPtheta1r2x1 
# "4/4"       "-16/4"        "16/4"        "16/4"       "-48/4"        "32/4" 
# vvPtheta2x2 vvPtheta2r1x1   vvPtheta2r2   vvPtheta3x1   vvPtheta3r1     vvPtheta4 
# "32/4"       "-64/4"        "16/4"        "32/4"       "-32/4"        "16/4" 
write_poly_formula(vv, round(res$coefficients * 4))
## simplified form
q_x_old <- q_x
q_x <- function(theta, r, x) {
  1/4 * (
  1*x^4 + -4*x^3*r + 4*r^2*x^2 + 4*x^3*theta +
    -12*x^2*theta*r + 8*r^2*theta*x + 8*theta^2*x^2 +
    -16*theta^2*r*x + 4*theta^2*r^2 + 8*theta^3*x +
    -8*theta^3*r + 4*theta^4)
}

## Simplify deriv_q_x
vv <- AlgDesign::gen.factorial(5, 3)
colnames(vv) <- c("theta", "r", "x")
V <- form_poly_matrix(vv, 4)
vals <- apply(vv, 1, function(v) deriv_q_x_raw(v[1], v[2], v[3]))
check_divis(vals, 64)
res <- lm(vals ~ 0+ V)
summary(res)
check_divis(res$coefficients, 1)
display_div(res$coefficients, 1)
# Vx3       Vr1x2       Vr2x1   Vtheta1x2 Vtheta1r1x1   Vtheta1r2   Vtheta2x1 
# "1/1"      "-3/1"       "2/1"       "3/1"      "-6/1"       "2/1"       "4/1" 
# Vtheta2r1     Vtheta3 
# "-4/1"       "2/1" 
deriv_q_x_coefs <- round(res$coefficients)
write_poly_formula(vv, deriv_q_x_coefs)
deriv_q_x <- function(theta, r, x) {
  1*x^3 + -3*x^2*r + 2*r^2*x + 3*x^2*theta + 
    -6*theta*r*x + 2*r^2*theta + 4*theta^2*x +
    -4*theta^2*r + 2*theta^3
}



## find coefs of second deriv

deriv2_q_x_raw <- function(theta, r, x)
  1e10 * (deriv_q_x(theta, r, x+1e-10) - deriv_q_x(theta, r, x))
vals <- apply(vv, 1, function(v) deriv2_q_x_raw(v[1], v[2], v[3]))
res <- lm(vals ~ 0+ V)
summary(res)
check_divis(res$coefficients, 1)
display_div(res$coefficients, 1)
# Vx2     Vr1x1       Vr2 Vtheta1x1 Vtheta1r1   Vtheta2 
# "3/1"    "-6/1"     "2/1"     "6/1"    "-6/1"     "4/1" 
deriv2_q_x_coefs <- round(res$coefficients)
write_poly_formula(vv, deriv2_q_x_coefs)
deriv2_q_x <- function(theta, r, x) {
  3*x^2 + -6*r*x + 2*r^2 + 6*theta*x + -6*theta*r + 4*theta^2
}

## check it out

theta <- runif(1); r <- runif(1); xs <- seq(0, 1, by = 1e-2)
x <- runif(1)
plot(xs, q_x(theta, r, xs), type = 'l')
lines(xs, q_x(theta, r, x) + (xs - x) * deriv_q_x(theta, r, x), col = "red")
lines(xs, q_x(theta, r, x) + (xs - x) * deriv_q_x(theta, r, x) +
        (xs - x)^2/2 * deriv2_q_x(theta, r, x), col = "blue")

#deriv_q_x(theta, r, x)
#numDeriv::grad(function(x) q_x(theta, r, x), x)
#deriv2_q_x(theta, r, x)
#numDeriv::hessian(function(x) q_x(theta, r, x), x)

####
##  Derivatives of ip(x)
####

deriv_sqrt_q <- function(theta, r, x)
  deriv_q_x(theta, r, x)/(2 * sqrt(q_x(theta, r, x)))
deriv2_sqrt_q <- function(theta, r, x) {
  qq <- q_x(theta, r, x)
  -1/4 * deriv_q_x(theta, r, x)^2 * qq^(-3/2) + 
    1/(2 * sqrt(qq)) * deriv2_q_x(theta, r, x)  
}

#deriv_sqrt_q(theta, r, x)
#numDeriv::grad(function(x) sqrt(q_x(theta, r, x)), x)
#deriv2_sqrt_q(theta, r, x)
#numDeriv::hessian(function(x) sqrt(q_x(theta, r, x)), x)

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

## check it out
theta <- 10 * runif(1); r <- theta * runif(1); xs <- seq(0, r, by = 1e-2)
x <- r * runif(1)
plot(xs, sqrt(q_x(theta, r, xs)), type = "l", main = "sqrt(q(x))")
lines(xs, sqrt(q_x(theta, r, x)) +
        (xs - x)*deriv_sqrt_q(theta, r, x), col = "red")
lines(xs, sqrt(q_x(theta, r, x)) +
        (xs - x)*deriv_sqrt_q(theta, r, x) +
        (xs - x)^2/2*deriv2_sqrt_q(theta, r, x), col = "blue")

theta <- 10 * runif(1); r <- theta * runif(1); xs <- seq(0, r, by = 1e-2)
x <- r * runif(1)
plot(xs, ip_x(theta, r, xs), type = "l", main = "ip(x)")
lines(xs, ip_x(theta, r, x) +
       (xs - x)*deriv_ip_x(theta, r, x), col = "red")
lines(xs, ip_x(theta, r, x) +
        (xs - x)*deriv_ip_x(theta, r, x) + 
        (xs - x)^2/2 * deriv2_ip_x(theta, r, x), col = "blue")

