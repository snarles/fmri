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
