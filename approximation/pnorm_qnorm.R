####
##  Formulas involving pnorm, qnorm
####

library(pracma); library(MASS)

## probability Pr[N(mu, sigma^2) < max_K N(0,1)]
normal_gumbel_formula <- function(mu, sigma2, k, reso = 1e4, naive = FALSE) {
  if (naive) {
    ind_lt <- rep(FALSE, reso)
    for (i in 1:reso) {
      zs <- rnorm(k)
      w <- rnorm(1, mean=mu, sd=sqrt(sigma2))
      ind_lt[i] <- (w < max(zs))
    }
    return(mean(ind_lt))
  }
  us <- (1:reso - .5)/reso
  mean(pnorm((qnorm(us^(1/k)) - mu)/sqrt(sigma2)))
}

## Probability that y_* < min_{i=1}^{K-1} y_i, where
## y_*    ( [a]  [b c c] )
## y_1 ~ N( [0], [c d e] )
## y_2    ( [0]  [c e d] )

pnormal_abcde <- function(aa, bb, cc, dd, ee, K, naive = TRUE) {
  if (naive) {
    p <- K
    mu <- rep(0, K); mu[1] <- aa
    Sigma <- diag(rep(dd - ee, K)) + ee
    Sigma[1, 1] <- bb
    Sigma[1, -1] <- cc; Sigma[-1, 1] <- cc
    stopifnot(min(eigen(Sigma)$values) > 0)
    stopifnot(dd > ee)
    stopifnot(ee > cc^2/bb)
    vals <- mvrnorm(1e6, mu, Sigma)
    pmins <- do.call(pmin, data.frame(vals))
    return(mean(vals[, 1] == pmins))
  }
  m <- -aa/sqrt(dd- ee)
  v <- (bb + ee - 2 * cc)/(dd - ee)
  1 - normal_gumbel_formula(m, v, K - 1)
}

####
##  TESTS
####

mu <- rnorm(1); sigma2 <- rexp(1); k <- 10
normal_gumbel_formula(mu, sigma2, k, naive = TRUE)
normal_gumbel_formula(mu, sigma2, k, naive = FALSE)

aa <- -3; bb <- 1.5; cc <- 2; dd <- 6; ee <- 3; K <- 5
pnormal_abcde(aa, bb, cc, dd, ee, K, TRUE)
pnormal_abcde(aa, bb, cc, dd, ee, K)
