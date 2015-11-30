####
##  Formulas involving pnorm, qnorm
####

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
## y_*
## y_1 ~ N(0, )
## y_2

####
##  TESTS
####

mu <- rnorm(1); sigma2 <- rexp(1); k <- 10
normal_gumbel_formula(mu, sigma2, k, naive = TRUE)
normal_gumbel_formula(mu, sigma2, k, naive = FALSE)
