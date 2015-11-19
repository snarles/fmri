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
  order1 <- eta - log(eta)/2 + (.25 * log(eta) - .5)/eta
  order2 <- order1 + (log(eta)^2/16 - 3/8 * log(eta) + 7/8)/(eta^2)
  order3 <- order2 + (log(eta)^3/48 - 7/32 * log(eta)^2 + 
                        (17/16) * log(eta) - 107/48)/(eta^3)
  list(true = true, o1 = sqrt(2 * order1), o2 = sqrt(2 * order2), 
       o3 = sqrt(2 * order3), eps = 1/eta)
}

qnorm2(0.9)
qnorm2(0.99)
qnorm2(0.999)
qnorm2(0.9999)
qnorm2(0.99999)



