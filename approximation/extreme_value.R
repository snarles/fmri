####
##  Does the min of ncx converge to a distribution?
####

## check the accuracy of the incomplete gamma for large shape params


k <- 1e5 # shape
q <- k  - 0.2 * sqrt(k)
pnorm(q=q, mean=k, sd=sqrt(k))  # should approach this as k to inf
pgamma(q=q, shape=k)


## find the median of a min_L chi^2(k, lambda)

## may not be necessary: pchisq is better!
ncx_cdf <- function(x, k, lambda, j.init = 0, j.end = 1000) {
  js <- j.init:j.end
  qs <- pgamma(x/2, (k + 2 * js)/2)
  weights <- dpois(js, lambda/2)
  sum(weights * qs)
}

pchisq(10, 5, 3)
ncx_cdf(10, 5, 3)

## log 1+x
log_1_plus <- function(x) {
  if (abs(x) > 1e-3) return(log(1 + x))
  sum(-1 * (-x)^(1:20)/(1:20))
}

med_loc <- function(k, L) qchisq(1 - (.5)^(1/L), k, k)

get_L <- function(x, k) {
  pp <- pchisq(x, k, k)
  log(1/2)/log_1_plus(-pp)
}

mincdf <- function(x, k, L) 1 - exp(L * log_1_plus(-pchisq(x, k, k)))
minpdf <- function(x, k, L) {
  pp <- pchisq(x, k, k)
  L * exp(L * log_1_plus(-pp))/(1 - pp)
}
  
  1 - exp(L * log_1_plus(-pchisq(x, k, k)))


L <- get_L(1, 10)
med_loc(10, L)

get_L(1, 1)


Ls <- matrix(0, 100, 10)
for (k in 1:10) {
  Ls[, k] <- sapply(1:100, function(i) get_L(k, i))  
}
matplot(log(Ls), type = "l")

Ls[100, 5]





mincdf(4, 100, Ls[100, 5])
mincdf(5, 100, Ls[100, 5])
mincdf(6, 100, Ls[100, 5])

minpdf(5, 80, Ls[80, 5])
minpdf(5, 100, Ls[100, 5])
