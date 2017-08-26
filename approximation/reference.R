####
##  Checking standard formulae
####

####
## Distribution of dot product
##   if t = <x, y> for x, y in S^{d-1}
##   then (t + 1)/2 ~ Beta((d-1)/2, (d-1)/2)
####

get_dot <- function(d) {
  x <- rnorm(d)
  x <- x/sqrt(sum(x^2))
  x[1]
}

mc.reps <- 1e4
d <- 10
ips <- sapply(1:mc.reps, function(i) get_dot(d))
plot(1:mc.reps/mc.reps, sort(pbeta((ips + 1)/2, (d-1)/2, (d-1)/2)), type = 'l')
abline(0, 1, col = "red")

####
## Chi distribution
####

library(Runuran)
dchi <- function(x, df) 2^(1 - (df/2)) * x^(df - 1) * exp(-x^2/2)/gamma(df/2)
x <- 20; df <- 10
ud(udchi(df), x)
dchi(x, df)

####
##  Normal quantile
####

qnorm2 <- function(p) {
  eta <- -log(2*sqrt(pi)*(1-p))
  true <- qnorm(p)
  order0 <- eta
  order0.5 <- order0 - log(eta)/2 
  order1 <- order0.5 + (.25 * log(eta) - .5)/eta
  order2 <- order1 + (log(eta)^2/16 - 3/8 * log(eta) + 7/8)/(eta^2)
  order3 <- order2 + (log(eta)^3/48 - 7/32 * log(eta)^2 + 
                        (17/16) * log(eta) - 107/48)/(eta^3)
  list(true = true, o0 = sqrt(2 * order0), o0.5 = sqrt(2 * order0.5),
       o1 = sqrt(2 * order1), o2 = sqrt(2 * order2), 
       o3 = sqrt(2 * order3), eps = 1/eta)
}

qnorm2(0.9)
qnorm2(0.99)
qnorm2(0.999)
qnorm2(0.9999)
qnorm2(0.99999)
qnorm2(0.999999)
qnorm2(0.9999999)


####
##  Exponential tilting for noncentral chi-squared
####

ncp <- rexp(1); df <- 10 * rexp(1); xs <- 1:5/10
tt <- -rexp(1)

quadfmla <- function(a, b, c, s = c(-1, 1))
  (-b + s*sqrt(b^2 - 4*a*c))/(2 * a)
psi <- function(tt) ncp * tt/(1 - 2*tt) - (df/2) * log(1-2*tt)
dpsi <- function(tt) (ncp + df)/(1 - 2*tt) + 2 * ncp * tt/(1-2*tt)^2
d2psi <- function(tt) 4*ncp/(1-2*tt)^2 + 8*ncp*tt/(1-2*tt)^3 + 2*df/(1-2*tt)^2
dpsi_inv <- function(x) quadfmla(4*x, -4*x+2*df, x - ncp - df, -1)

exp(tt * xs - psi(tt)) * dchisq(xs, df, ncp)
(1-2*tt) * dchisq(xs * (1-2*tt), df, ncp/(1-2*tt))

x <- dpsi(0) * 0.1
tt <- dpsi_inv(x)

pchisq(x, df, ncp)

reso = 1e4
sig <- 1/(1-2*tt)
ncp2 <- ncp*sig

xs <- x * (1:reso)/reso
dlta <- xs[2] - xs[1]
exp(psi(tt)) * dlta * 
  sum(exp(-tt * xs) * (1/sig) * dchisq(xs/sig, df, ncp2))

plot(xs, exp(-tt * xs) * (1/sig) * dchisq(xs/sig, df, ncp2), type = "l")
plot(xs, dchisq(xs/sig, df, ncp2), type = "l")
