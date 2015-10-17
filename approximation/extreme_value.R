####
##  Does the min of ncx converge to a distribution?
####

## check the accuracy of the incomplete gamma for large shape params


k <- 1e5 # shape
q <- k  - 0.2 * sqrt(k)
pnorm(q=q, mean=k, sd=sqrt(k))  # should approach this as k to inf
pgamma(q=q, shape=k)


## find the median of a min_L chi^2(k, lambda)

ncx_cdf <- function(x, k, lambda, j.init = 0, j.end = 1000) {
  js <- j.init:j.end
  qs <- pgamma(x/2, (k + 2 * js)/2)
  weights <- dpois(js, lambda/2)
  sum(weights * qs)
}

pchisq(10, 5, 3)
ncx_cdf(10, 5, 3)


