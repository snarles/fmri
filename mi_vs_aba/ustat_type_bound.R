## Get a bound on expected max vs var max using ustatistic type thing

compute_moments <- function(as, bs) {
  k <- length(as) - 1
  ss <- sum(as)
  ms <- ss - as + bs
  c(mean(ms), mean(ms^2))
}



compute_moments(c(0.6, 0.7, 0.8), c(0.2, 0.2, 0.3))

## explore empirical dists without optimization

k <- 4
pts <- do.call(rbind, lapply(1:1e5, function(i) {
  as <- runif(k + 1)
  bs <- runif(k + 1) * as
  sb <- sum(bs)
  sa <- sum(as)
  if (sb > (k + 1) - sa) {
    return(c(NA, NA))
    ##bs <- bs * ((k + 1) - sa)/sb
  } else if (sb < 1 - (sa)/(k+1)) {
    return(c(NA, NA))
  }
  compute_moments(as, bs)
}))
plot(pts, pch = ".")
lines(seq(0, k, 0.1), seq(0, k, 0.1)^2, col = "red")
