### Fix v_1 < ... < v_d
### Let X_i = v_{pi_i} (single X)
### Get coordinatewise max of k independent X's, mean and var of L1-norm/k??


library(combinat)
library(AlgDesign)

compute_weights <- function(d, k, mc.reps = 1e4) {
  if (mc.reps == Inf) {
    perms <- do.call(rbind, permn(1:d))
    np <- nrow(perms)
    kc <- gen.factorial(rep(np, k), center = FALSE)
    bigmat <- array(t(perms)[, as.numeric(t(kc))], dim = c(d, k, nrow(kc)))
    lala <- apply(bigmat, c(1, 3), max)
  } else {
    lala <- sapply(1:mc.reps, function(i) {
      pis <- sapply(1:k, function(i) sample(d, d))
      apply(pis, 1, max)
    })
  }
  ncts <- sapply(1:d, function(i) colSums(lala == i))
  wts <- colSums(ncts)
  wts <- wts/sum(wts)
  ## now build quadratic stats
  wmat <- matrix(0, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      wmat[i, j] <- sum(ncts[, i] * ncts[, j])
    }
  }
  wmat <- wmat/sum(wmat)
  wmat <- wmat - wts %*% t(wts)
  return(list(wts = wts, wmat = wmat))
}

compute_moments <- function(v, wtz) {
  v <- sort(v, decreasing = FALSE)
  mu <- sum(v * wtz$wts)
  mu2 <- t(v) %*% wtz$wmat %*% v
  c(mu = mu, sig = sqrt(mu2))
}

compute_moments_naive <- function(v, k, mc.reps = 1e4) {
  d <- length(v)
  lala <- sapply(1:mc.reps, function(i) {
    pis <- sapply(1:k, function(i) sample(v, d))
    apply(pis, 1, max)
  })
  lal2 <- colMeans(lala)
  c(mu = mean(lal2), sig = sd(lal2))
}

# compute_weights(4, 3, mc.reps = 1e2)
# compute_weights(4, 3, mc.reps = 1e4)
# compute_weights(4, 3, mc.reps = Inf)
# 
# v <- runif(3)
# v <- v/sum(v)
# 
# k <- 3
# wtz <- compute_weights(d = length(v), k, mc.reps = 1e4)
# 
# compute_moments_naive(v, k)
# compute_moments(v, wtz)

